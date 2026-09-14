#include "webp_backend.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <queue>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::graphics {
namespace {

class BitWriter final {
public:
    void writeBits(std::uint32_t value, unsigned count) {
        bitBuffer_ |= static_cast<std::uint64_t>(value) << bitCount_;
        bitCount_ += count;
        while (bitCount_ >= 8u) {
            bytes_.push_back(static_cast<std::uint8_t>(bitBuffer_ & 0xffu));
            bitBuffer_ >>= 8u;
            bitCount_ -= 8u;
        }
    }

    [[nodiscard]] std::vector<std::uint8_t> finish() {
        if (bitCount_ != 0u)
            bytes_.push_back(static_cast<std::uint8_t>(bitBuffer_ & 0xffu));
        bitBuffer_ = 0;
        bitCount_ = 0;
        return std::move(bytes_);
    }

private:
    std::vector<std::uint8_t> bytes_;
    std::uint64_t bitBuffer_ = 0;
    unsigned bitCount_ = 0;
};

[[nodiscard]] std::uint32_t reverseBits(std::uint32_t value, unsigned count) noexcept {
    std::uint32_t reversed = 0;
    for (unsigned i = 0; i < count; ++i) {
        reversed = (reversed << 1u) | (value & 1u);
        value >>= 1u;
    }
    return reversed;
}

void appendLittleEndian32(std::string& output, std::uint32_t value) {
    output.push_back(static_cast<char>(value & 0xffu));
    output.push_back(static_cast<char>((value >> 8u) & 0xffu));
    output.push_back(static_cast<char>((value >> 16u) & 0xffu));
    output.push_back(static_cast<char>((value >> 24u) & 0xffu));
}

struct PrefixCode final {
    std::uint16_t code = 0;
    std::uint8_t bits = 0;
};

// 頻度からcanonical prefix codeのbit-lengthを作る。
// VP8L data treeは最大15 bit，code-length treeは最大7 bit。
// 通常のHuffman treeが上限を越える場合だけ，完全木を保つbalanced fallbackへ落とす。
[[nodiscard]] std::vector<std::uint8_t> makeHuffmanLengths(
    const std::vector<std::uint64_t>& frequencies,
    unsigned maxBits) {
    std::vector<std::uint8_t> lengths(frequencies.size(), static_cast<std::uint8_t>(0));
    if (frequencies.empty() || maxBits == 0u)
        return lengths;

    struct ActiveSymbol final {
        std::uint64_t frequency = 0;
        unsigned symbol = 0;
    };
    std::vector<ActiveSymbol> active;
    active.reserve(frequencies.size());
    for (unsigned symbol = 0; symbol < frequencies.size(); ++symbol) {
        if (frequencies[symbol] != 0u)
            active.push_back(ActiveSymbol{frequencies[symbol], symbol});
    }

    // VP8Lではempty treeもsingle leafとしてsymbol 0を記述できる。
    if (active.empty()) {
        if (!lengths.empty())
            lengths[0] = 1u;
        return lengths;
    }
    if (active.size() == 1u) {
        lengths[active.front().symbol] = 1u;
        return lengths;
    }

    struct Node final {
        std::uint64_t frequency = 0;
        unsigned minSymbol = 0;
        int left = -1;
        int right = -1;
        int symbol = -1;
    };
    struct QueueItem final {
        std::uint64_t frequency = 0;
        unsigned minSymbol = 0;
        int node = -1;
    };
    struct QueueGreater final {
        [[nodiscard]] bool operator()(const QueueItem& a, const QueueItem& b) const noexcept {
            if (a.frequency != b.frequency)
                return a.frequency > b.frequency;
            if (a.minSymbol != b.minSymbol)
                return a.minSymbol > b.minSymbol;
            return a.node > b.node;
        }
    };

    std::vector<Node> nodes;
    nodes.reserve(active.size() * 2u - 1u);
    std::priority_queue<QueueItem, std::vector<QueueItem>, QueueGreater> queue;
    for (const ActiveSymbol item : active) {
        const int index = static_cast<int>(nodes.size());
        nodes.push_back(Node{item.frequency, item.symbol, -1, -1,
            static_cast<int>(item.symbol)});
        queue.push(QueueItem{item.frequency, item.symbol, index});
    }

    while (queue.size() > 1u) {
        const QueueItem a = queue.top();
        queue.pop();
        const QueueItem b = queue.top();
        queue.pop();
        const int index = static_cast<int>(nodes.size());
        nodes.push_back(Node{
            a.frequency + b.frequency,
            std::min(a.minSymbol, b.minSymbol),
            a.node,
            b.node,
            -1});
        queue.push(QueueItem{
            a.frequency + b.frequency,
            std::min(a.minSymbol, b.minSymbol),
            index});
    }

    unsigned deepest = 0u;
    std::vector<std::pair<int, unsigned>> stack;
    stack.emplace_back(queue.top().node, 0u);
    while (!stack.empty()) {
        const auto [nodeIndex, depth] = stack.back();
        stack.pop_back();
        const Node& node = nodes[static_cast<std::size_t>(nodeIndex)];
        if (node.symbol >= 0) {
            const unsigned leafDepth = std::max(1u, depth);
            lengths[static_cast<unsigned>(node.symbol)] =
                static_cast<std::uint8_t>(leafDepth);
            deepest = std::max(deepest, leafDepth);
            continue;
        }
        stack.emplace_back(node.left, depth + 1u);
        stack.emplace_back(node.right, depth + 1u);
    }
    if (deepest <= maxBits)
        return lengths;

    // 極端に偏った頻度で制限を越えた場合の保険。
    // n leafの完全木をfloor/ceil(log2 n)で構成し，頻出symbolへ短いcodeを割り当てる。
    std::fill(lengths.begin(), lengths.end(), static_cast<std::uint8_t>(0));
    std::sort(active.begin(), active.end(), [](const ActiveSymbol& a, const ActiveSymbol& b) {
        if (a.frequency != b.frequency)
            return a.frequency > b.frequency;
        return a.symbol < b.symbol;
    });
    unsigned longBits = 1u;
    while ((std::size_t{1} << longBits) < active.size())
        ++longBits;
    if (longBits > maxBits)
        return {};
    const unsigned shortBits = longBits - 1u;
    const std::size_t shortCount = (std::size_t{1} << longBits) - active.size();
    for (std::size_t i = 0; i < active.size(); ++i) {
        lengths[active[i].symbol] = static_cast<std::uint8_t>(
            i < shortCount ? shortBits : longBits);
    }
    return lengths;
}

[[nodiscard]] std::vector<PrefixCode> makeCanonicalCodes(
    const std::vector<std::uint8_t>& lengths) {
    std::array<unsigned, 16> count{};
    unsigned nonzero = 0u;
    unsigned singleSymbol = 0u;
    for (unsigned symbol = 0; symbol < lengths.size(); ++symbol) {
        const std::uint8_t length = lengths[symbol];
        if (length > 15u)
            return {};
        if (length != 0u) {
            ++count[length];
            ++nonzero;
            singleSymbol = symbol;
        }
    }

    std::vector<PrefixCode> table(lengths.size());
    if (nonzero == 1u) {
        // single-leaf treeはcode length 1として伝送するが，data bitは消費しない。
        table[singleSymbol] = PrefixCode{0u, 0u};
        return table;
    }

    std::array<unsigned, 16> next{};
    unsigned code = 0;
    for (unsigned bits = 1; bits <= 15; ++bits) {
        code = (code + count[bits - 1]) << 1u;
        next[bits] = code;
    }

    for (std::size_t symbol = 0; symbol < lengths.size(); ++symbol) {
        const unsigned bits = lengths[symbol];
        if (bits == 0u)
            continue;
        table[symbol] = PrefixCode{
            static_cast<std::uint16_t>(next[bits]++),
            static_cast<std::uint8_t>(bits)};
    }
    return table;
}

void writePrefixSymbol(
    BitWriter& writer,
    const std::vector<PrefixCode>& table,
    unsigned symbol) {
    const PrefixCode code = table[symbol];
    writer.writeBits(reverseBits(code.code, code.bits), code.bits);
}

struct CodeLengthToken final {
    unsigned symbol = 0u;
    unsigned extraBits = 0u;
    unsigned extraValue = 0u;
};

[[nodiscard]] std::vector<CodeLengthToken> makeCodeLengthTokens(
    const std::vector<std::uint8_t>& lengths) {
    std::vector<CodeLengthToken> tokens;
    tokens.reserve(lengths.size());

    std::size_t position = 0u;
    while (position < lengths.size()) {
        const unsigned value = lengths[position];
        std::size_t run = 1u;
        while (position + run < lengths.size() && lengths[position + run] == value)
            ++run;

        if (value == 0u) {
            std::size_t remaining = run;
            while (remaining >= 11u) {
                std::size_t count = std::min<std::size_t>(remaining, 138u);
                const std::size_t tail = remaining - count;
                if (tail == 1u && count >= 13u)
                    count -= 2u;
                else if (tail == 2u && count >= 12u)
                    count -= 1u;
                tokens.push_back(CodeLengthToken{18u, 7u,
                    static_cast<unsigned>(count - 11u)});
                remaining -= count;
            }
            while (remaining >= 3u) {
                std::size_t count = std::min<std::size_t>(remaining, 10u);
                const std::size_t tail = remaining - count;
                if (tail == 1u && count >= 5u)
                    count -= 2u;
                else if (tail == 2u && count >= 4u)
                    count -= 1u;
                tokens.push_back(CodeLengthToken{17u, 3u,
                    static_cast<unsigned>(count - 3u)});
                remaining -= count;
            }
            while (remaining-- != 0u)
                tokens.push_back(CodeLengthToken{0u, 0u, 0u});
        }
        else {
            tokens.push_back(CodeLengthToken{value, 0u, 0u});
            std::size_t remaining = run - 1u;
            while (remaining >= 3u) {
                std::size_t count = std::min<std::size_t>(remaining, 6u);
                const std::size_t tail = remaining - count;
                if (tail == 1u && count >= 5u)
                    count -= 2u;
                else if (tail == 2u && count >= 4u)
                    count -= 1u;
                tokens.push_back(CodeLengthToken{16u, 2u,
                    static_cast<unsigned>(count - 3u)});
                remaining -= count;
            }
            while (remaining-- != 0u)
                tokens.push_back(CodeLengthToken{value, 0u, 0u});
        }
        position += run;
    }
    return tokens;
}

void writePrefixCodeLengths(BitWriter& writer, const std::vector<std::uint8_t>& lengths) {
    std::vector<unsigned> activeSymbols;
    for (unsigned symbol = 0; symbol < lengths.size(); ++symbol) {
        if (lengths[symbol] != 0u)
            activeSymbols.push_back(symbol);
    }

    // 1～2 leafかつsymbol <=255ならsimple formが最小で，single leafはdata bit不要。
    if (!activeSymbols.empty() && activeSymbols.size() <= 2u
        && activeSymbols.back() <= 255u) {
        writer.writeBits(1u, 1); // simple code length code
        writer.writeBits(static_cast<unsigned>(activeSymbols.size() - 1u), 1);
        const unsigned first = activeSymbols[0];
        if (first <= 1u) {
            writer.writeBits(0u, 1); // first symbol uses 1 bit
            writer.writeBits(first, 1);
        }
        else {
            writer.writeBits(1u, 1); // first symbol uses 8 bits
            writer.writeBits(first, 8);
        }
        if (activeSymbols.size() == 2u)
            writer.writeBits(activeSymbols[1], 8);
        return;
    }

    const std::vector<CodeLengthToken> tokens = makeCodeLengthTokens(lengths);
    std::vector<std::uint64_t> tokenFrequencies(19u, 0u);
    for (const CodeLengthToken& token : tokens)
        ++tokenFrequencies[token.symbol];
    const std::vector<std::uint8_t> codeLengthLengths =
        makeHuffmanLengths(tokenFrequencies, 7u);
    const auto codeLengthCodes = makeCanonicalCodes(codeLengthLengths);

    writer.writeBits(0u, 1); // normal code length code
    static constexpr std::array<unsigned, 19> order{
        17,18,0,1,2,3,4,5,16,6,7,8,9,10,11,12,13,14,15};
    unsigned numCodeLengths = 4u;
    for (unsigned i = 0; i < order.size(); ++i) {
        if (codeLengthLengths[order[i]] != 0u)
            numCodeLengths = std::max(numCodeLengths, i + 1u);
    }
    writer.writeBits(numCodeLengths - 4u, 4);
    for (unsigned i = 0; i < numCodeLengths; ++i)
        writer.writeBits(codeLengthLengths[order[i]], 3);

    // alphabet全体を記述する。末尾zeroは17/18で圧縮するためcustom max_symbolは不要。
    writer.writeBits(0u, 1);
    for (const CodeLengthToken& token : tokens) {
        writePrefixSymbol(writer, codeLengthCodes, token.symbol);
        if (token.extraBits != 0u)
            writer.writeBits(token.extraValue, token.extraBits);
    }
}

struct PrefixInteger final {
    unsigned code = 0;
    unsigned extraBits = 0;
    unsigned extraValue = 0;
};

[[nodiscard]] PrefixInteger encodePrefixInteger(unsigned value) noexcept {
    if (value <= 4u)
        return PrefixInteger{value - 1u, 0u, 0u};

    for (unsigned code = 4; code < 40; ++code) {
        const unsigned extraBits = (code - 2u) >> 1u;
        const unsigned offset = (2u + (code & 1u)) << extraBits;
        const unsigned first = offset + 1u;
        const unsigned last = offset + (1u << extraBits);
        if (value >= first && value <= last)
            return PrefixInteger{code, extraBits, value - first};
    }
    return PrefixInteger{39u, 18u, 0u};
}

[[nodiscard]] std::uint32_t pixelAt(
    const std::vector<std::uint8_t>& rgba,
    std::size_t position) noexcept {
    const std::size_t i = position * 4u;
    return static_cast<std::uint32_t>(rgba[i + 0])
        | (static_cast<std::uint32_t>(rgba[i + 1]) << 8u)
        | (static_cast<std::uint32_t>(rgba[i + 2]) << 16u)
        | (static_cast<std::uint32_t>(rgba[i + 3]) << 24u);
}

// VP8Lのcolor cache hashは仕様上のARGB 32-bit値
// (A:31..24, R:23..16, G:15..8, B:7..0) に対して行う。
[[nodiscard]] std::uint32_t webpArgbAt(
    const std::vector<std::uint8_t>& rgba,
    std::size_t position) noexcept {
    const std::size_t i = position * 4u;
    return static_cast<std::uint32_t>(rgba[i + 2u])
        | (static_cast<std::uint32_t>(rgba[i + 1u]) << 8u)
        | (static_cast<std::uint32_t>(rgba[i + 0u]) << 16u)
        | (static_cast<std::uint32_t>(rgba[i + 3u]) << 24u);
}

class ColorCache final {
public:
    explicit ColorCache(unsigned bits)
        : bits_(bits), entries_(std::size_t{1} << bits, 0u) {}

    [[nodiscard]] unsigned index(std::uint32_t argb) const noexcept {
        constexpr std::uint32_t multiplier = 0x1e35a7bdu;
        return static_cast<unsigned>((argb * multiplier) >> (32u - bits_));
    }

    [[nodiscard]] std::optional<unsigned> find(std::uint32_t argb) const noexcept {
        const unsigned key = index(argb);
        if (entries_[key] == argb)
            return key;
        return std::nullopt;
    }

    void insert(std::uint32_t argb) noexcept {
        entries_[index(argb)] = argb;
    }

private:
    unsigned bits_ = 0;
    std::vector<std::uint32_t> entries_;
};

[[nodiscard]] std::vector<std::uint8_t> makeProvisionalGreenLengths(unsigned colorCacheSize) {
    const unsigned alphabetSize = 280u + colorCacheSize;
    unsigned longBits = 1u;
    while ((1u << longBits) < alphabetSize)
        ++longBits;
    const unsigned shortBits = longBits - 1u;
    const unsigned shortCount = (1u << longBits) - alphabetSize;

    std::vector<std::uint8_t> lengths(alphabetSize, static_cast<std::uint8_t>(longBits));
    unsigned remaining = shortCount;

    // backward-reference length codeをまず短くする。次にcolor cache codeを優先する。
    // 残りがあれば低green literalへ割り当てる。完全二分木のbit-length集合は維持される。
    for (unsigned symbol = 256u; symbol < 280u && remaining != 0u; ++symbol, --remaining)
        lengths[symbol] = static_cast<std::uint8_t>(shortBits);
    for (unsigned symbol = 280u; symbol < alphabetSize && remaining != 0u; ++symbol, --remaining)
        lengths[symbol] = static_cast<std::uint8_t>(shortBits);
    for (unsigned symbol = 0u; symbol < 256u && remaining != 0u; ++symbol, --remaining)
        lengths[symbol] = static_cast<std::uint8_t>(shortBits);
    return lengths;
}

[[nodiscard]] std::uint16_t hash3Pixels(
    const std::vector<std::uint8_t>& rgba,
    std::size_t position) noexcept {
    std::uint32_t value = pixelAt(rgba, position) * 0x1e35a7bdu;
    value ^= pixelAt(rgba, position + 1u) * 0x9e3779b1u;
    value ^= pixelAt(rgba, position + 2u) * 0x85ebca6bu;
    value ^= value >> 16u;
    return static_cast<std::uint16_t>(value & 0xffffu);
}

struct Match final {
    std::size_t length = 0;
    std::size_t distance = 0;
};

struct Pixel final {
    std::uint8_t red = 0;
    std::uint8_t green = 0;
    std::uint8_t blue = 0;
    std::uint8_t alpha = 0;
};

[[nodiscard]] Pixel pixelAtRgba(
    const std::vector<std::uint8_t>& rgba,
    unsigned width,
    unsigned x,
    unsigned y) noexcept {
    const std::size_t i =
        (static_cast<std::size_t>(y) * width + x) * 4u;
    return Pixel{
        rgba[i + 0u], rgba[i + 1u], rgba[i + 2u], rgba[i + 3u]};
}

[[nodiscard]] std::uint8_t average2(std::uint8_t a, std::uint8_t b) noexcept {
    return static_cast<std::uint8_t>(
        (static_cast<unsigned>(a) + static_cast<unsigned>(b)) >> 1u);
}

[[nodiscard]] Pixel average2(Pixel a, Pixel b) noexcept {
    return Pixel{
        average2(a.red, b.red),
        average2(a.green, b.green),
        average2(a.blue, b.blue),
        average2(a.alpha, b.alpha)};
}

[[nodiscard]] std::uint8_t clampByte(int value) noexcept {
    if (value < 0)
        return 0u;
    if (value > 255)
        return 255u;
    return static_cast<std::uint8_t>(value);
}

[[nodiscard]] Pixel clampAddSubtractFull(Pixel left, Pixel top, Pixel topLeft) noexcept {
    return Pixel{
        clampByte(static_cast<int>(left.red) + top.red - topLeft.red),
        clampByte(static_cast<int>(left.green) + top.green - topLeft.green),
        clampByte(static_cast<int>(left.blue) + top.blue - topLeft.blue),
        clampByte(static_cast<int>(left.alpha) + top.alpha - topLeft.alpha)};
}

[[nodiscard]] std::uint8_t clampAddSubtractHalf(
    std::uint8_t a,
    std::uint8_t b) noexcept {
    const int ai = static_cast<int>(a);
    return clampByte(ai + (ai - static_cast<int>(b)) / 2);
}

[[nodiscard]] Pixel clampAddSubtractHalf(Pixel a, Pixel b) noexcept {
    return Pixel{
        clampAddSubtractHalf(a.red, b.red),
        clampAddSubtractHalf(a.green, b.green),
        clampAddSubtractHalf(a.blue, b.blue),
        clampAddSubtractHalf(a.alpha, b.alpha)};
}

[[nodiscard]] int absInt(int value) noexcept {
    return value < 0 ? -value : value;
}

// VP8L predictor mode 11。componentごとではなく，ARGB全体のManhattan距離で
// left/topのどちらか一方のpixelを選ぶ。
[[nodiscard]] Pixel selectPredictor(Pixel left, Pixel top, Pixel topLeft) noexcept {
    const int estimateAlpha = static_cast<int>(left.alpha) + top.alpha - topLeft.alpha;
    const int estimateRed = static_cast<int>(left.red) + top.red - topLeft.red;
    const int estimateGreen = static_cast<int>(left.green) + top.green - topLeft.green;
    const int estimateBlue = static_cast<int>(left.blue) + top.blue - topLeft.blue;

    const int leftDistance =
        absInt(estimateAlpha - left.alpha)
        + absInt(estimateRed - left.red)
        + absInt(estimateGreen - left.green)
        + absInt(estimateBlue - left.blue);
    const int topDistance =
        absInt(estimateAlpha - top.alpha)
        + absInt(estimateRed - top.red)
        + absInt(estimateGreen - top.green)
        + absInt(estimateBlue - top.blue);
    return leftDistance < topDistance ? left : top;
}

[[nodiscard]] Pixel predictorForMode(
    unsigned mode,
    Pixel left,
    Pixel top,
    Pixel topRight,
    Pixel topLeft) noexcept {
    switch (mode) {
    case 0:
        return Pixel{0u, 0u, 0u, 255u};
    case 1:
        return left;
    case 2:
        return top;
    case 3:
        return topRight;
    case 4:
        return topLeft;
    case 5:
        return average2(average2(left, topRight), top);
    case 6:
        return average2(left, topLeft);
    case 7:
        return average2(left, top);
    case 8:
        return average2(topLeft, top);
    case 9:
        return average2(top, topRight);
    case 10:
        return average2(average2(left, topLeft), average2(top, topRight));
    case 11:
        return selectPredictor(left, top, topLeft);
    case 12:
        return clampAddSubtractFull(left, top, topLeft);
    case 13:
        return clampAddSubtractHalf(average2(left, top), topLeft);
    default:
        return Pixel{};
    }
}

[[nodiscard]] Pixel predictorAt(
    const std::vector<std::uint8_t>& rgba,
    unsigned width,
    unsigned x,
    unsigned y,
    unsigned mode) noexcept {
    if (x == 0u && y == 0u)
        return Pixel{0u, 0u, 0u, 255u};
    if (y == 0u)
        return pixelAtRgba(rgba, width, x - 1u, y);
    if (x == 0u)
        return pixelAtRgba(rgba, width, x, y - 1u);

    const Pixel left = pixelAtRgba(rgba, width, x - 1u, y);
    const Pixel top = pixelAtRgba(rgba, width, x, y - 1u);
    const Pixel topLeft = pixelAtRgba(rgba, width, x - 1u, y - 1u);
    const Pixel topRight = x + 1u < width
        ? pixelAtRgba(rgba, width, x + 1u, y - 1u)
        // VP8L仕様では右端だけ，同じ行の左端pixelをTRとして使う。
        : pixelAtRgba(rgba, width, 0u, y);
    return predictorForMode(mode, left, top, topRight, topLeft);
}

[[nodiscard]] std::uint8_t residualByte(
    std::uint8_t actual,
    std::uint8_t predicted) noexcept {
    return static_cast<std::uint8_t>(
        static_cast<unsigned>(actual) - static_cast<unsigned>(predicted));
}

[[nodiscard]] unsigned residualCostByte(std::uint8_t residual) noexcept {
    if (residual == 0u)
        return 0u;
    const unsigned value = residual;
    const unsigned distanceToZero = std::min(value, 256u - value);
    // exact zeroを強く優先しつつ，小さいmod-256残差を次点にする。
    return 4u + distanceToZero;
}

[[nodiscard]] std::uint64_t predictorModeCost(
    const std::vector<std::uint8_t>& rgba,
    unsigned width,
    unsigned xBegin,
    unsigned yBegin,
    unsigned xEnd,
    unsigned yEnd,
    unsigned mode) noexcept {
    std::uint64_t cost = 0;
    for (unsigned y = yBegin; y < yEnd; ++y) {
        for (unsigned x = xBegin; x < xEnd; ++x) {
            // top row / left columnはtransform modeに依存しないため選択costから除く。
            if (x == 0u || y == 0u)
                continue;
            const Pixel actual = pixelAtRgba(rgba, width, x, y);
            const Pixel predicted = predictorAt(rgba, width, x, y, mode);
            cost += residualCostByte(residualByte(actual.red, predicted.red));
            cost += residualCostByte(residualByte(actual.green, predicted.green));
            cost += residualCostByte(residualByte(actual.blue, predicted.blue));
            cost += residualCostByte(residualByte(actual.alpha, predicted.alpha));
        }
    }
    return cost;
}

struct PredictorTransform final {
    unsigned sizeBits = 4u;
    unsigned transformWidth = 0u;
    unsigned transformHeight = 0u;
    std::vector<std::uint8_t> modesRgba;
    std::vector<std::uint8_t> residualRgba;
};

[[nodiscard]] PredictorTransform makePredictorTransform(
    const RasterImage& image,
    unsigned sizeBits) {
    PredictorTransform transform;
    transform.sizeBits = sizeBits;
    const unsigned blockSize = 1u << sizeBits;
    transform.transformWidth = (image.widthPx + blockSize - 1u) / blockSize;
    transform.transformHeight = (image.heightPx + blockSize - 1u) / blockSize;
    transform.modesRgba.resize(
        static_cast<std::size_t>(transform.transformWidth) * transform.transformHeight * 4u);
    transform.residualRgba.resize(image.rgba.size());

    std::vector<std::uint8_t> modes(
        static_cast<std::size_t>(transform.transformWidth) * transform.transformHeight, 1u);

    // 16x16 blockごとに14 predictorを全探索する。mode map自体は小さく，
    // transform選択をdeterministicかつboundedに保てる。
    for (unsigned by = 0; by < transform.transformHeight; ++by) {
        const unsigned yBegin = by * blockSize;
        const unsigned yEnd = std::min(image.heightPx, yBegin + blockSize);
        for (unsigned bx = 0; bx < transform.transformWidth; ++bx) {
            const unsigned xBegin = bx * blockSize;
            const unsigned xEnd = std::min(image.widthPx, xBegin + blockSize);

            unsigned bestMode = 0u;
            std::uint64_t bestCost = std::numeric_limits<std::uint64_t>::max();
            for (unsigned mode = 0; mode < 14u; ++mode) {
                const std::uint64_t cost = predictorModeCost(
                    image.rgba,
                    image.widthPx,
                    xBegin,
                    yBegin,
                    xEnd,
                    yEnd,
                    mode);
                if (cost < bestCost) {
                    bestCost = cost;
                    bestMode = mode;
                }
            }
            modes[static_cast<std::size_t>(by) * transform.transformWidth + bx] =
                static_cast<std::uint8_t>(bestMode);
        }
    }

    for (std::size_t i = 0; i < modes.size(); ++i) {
        const std::size_t p = i * 4u;
        transform.modesRgba[p + 0u] = 0u;
        transform.modesRgba[p + 1u] = modes[i];
        transform.modesRgba[p + 2u] = 0u;
        transform.modesRgba[p + 3u] = 255u;
    }

    for (unsigned y = 0; y < image.heightPx; ++y) {
        for (unsigned x = 0; x < image.widthPx; ++x) {
            const unsigned bx = x >> sizeBits;
            const unsigned by = y >> sizeBits;
            const unsigned mode = modes[
                static_cast<std::size_t>(by) * transform.transformWidth + bx];
            const Pixel actual = pixelAtRgba(image.rgba, image.widthPx, x, y);
            const Pixel predicted = predictorAt(
                image.rgba, image.widthPx, x, y, mode);
            const std::size_t p =
                (static_cast<std::size_t>(y) * image.widthPx + x) * 4u;
            transform.residualRgba[p + 0u] = residualByte(actual.red, predicted.red);
            transform.residualRgba[p + 1u] = residualByte(actual.green, predicted.green);
            transform.residualRgba[p + 2u] = residualByte(actual.blue, predicted.blue);
            transform.residualRgba[p + 3u] = residualByte(actual.alpha, predicted.alpha);
        }
    }

    return transform;
}

// VP8Lのdistance code 1..120は，scan-line距離そのものではなく
// 現在pixel近傍の2D offsetを短いcodeへ写像する。仕様の順序をそのまま保持する。
struct DistanceOffset final {
    int x = 0;
    int y = 0;
};

static constexpr std::array<DistanceOffset, 120> kDistanceMap{{
    { 0, 1}, { 1, 0}, { 1, 1}, {-1, 1}, { 0, 2}, { 2, 0}, { 1, 2}, {-1, 2},
    { 2, 1}, {-2, 1}, { 2, 2}, {-2, 2}, { 0, 3}, { 3, 0}, { 1, 3}, {-1, 3},
    { 3, 1}, {-3, 1}, { 2, 3}, {-2, 3}, { 3, 2}, {-3, 2}, { 0, 4}, { 4, 0},
    { 1, 4}, {-1, 4}, { 4, 1}, {-4, 1}, { 3, 3}, {-3, 3}, { 2, 4}, {-2, 4},
    { 4, 2}, {-4, 2}, { 0, 5}, { 3, 4}, {-3, 4}, { 4, 3}, {-4, 3}, { 5, 0},
    { 1, 5}, {-1, 5}, { 5, 1}, {-5, 1}, { 2, 5}, {-2, 5}, { 5, 2}, {-5, 2},
    { 4, 4}, {-4, 4}, { 3, 5}, {-3, 5}, { 5, 3}, {-5, 3}, { 0, 6}, { 6, 0},
    { 1, 6}, {-1, 6}, { 6, 1}, {-6, 1}, { 2, 6}, {-2, 6}, { 6, 2}, {-6, 2},
    { 4, 5}, {-4, 5}, { 5, 4}, {-5, 4}, { 3, 6}, {-3, 6}, { 6, 3}, {-6, 3},
    { 0, 7}, { 7, 0}, { 1, 7}, {-1, 7}, { 5, 5}, {-5, 5}, { 7, 1}, {-7, 1},
    { 4, 6}, {-4, 6}, { 6, 4}, {-6, 4}, { 2, 7}, {-2, 7}, { 7, 2}, {-7, 2},
    { 3, 7}, {-3, 7}, { 7, 3}, {-7, 3}, { 5, 6}, {-5, 6}, { 6, 5}, {-6, 5},
    { 8, 0}, { 4, 7}, {-4, 7}, { 7, 4}, {-7, 4}, { 8, 1}, { 8, 2}, { 6, 6},
    {-6, 6}, { 8, 3}, { 5, 7}, {-5, 7}, { 7, 5}, {-7, 5}, { 8, 4}, { 6, 7},
    {-6, 7}, { 7, 6}, {-7, 6}, { 8, 5}, { 7, 7}, {-7, 7}, { 8, 6}, { 8, 7}
}};

[[nodiscard]] unsigned distanceBitCost(
    unsigned mappedDistance,
    const std::vector<PrefixCode>& distanceCodes) noexcept {
    const PrefixInteger prefix = encodePrefixInteger(mappedDistance);
    return distanceCodes[prefix.code].bits + prefix.extraBits;
}

// LZ77 match finderはscan-line pixel distanceを返す。
// VP8Lのspecial 2D codeへ逆写像し，現在のdistance prefix treeで最小bit costとなる
// codeを選ぶ。幅が小さい画像では複数の2D offsetが同じscan distanceへ潰れるため，
// 最初のcodeを盲目的に採用せずcost比較する。
[[nodiscard]] std::vector<unsigned> makeMappedDistanceTable(
    unsigned imageWidth,
    const std::vector<PrefixCode>& distanceCodes) {
    std::vector<unsigned> table(65537u, 0u);
    for (unsigned distance = 1; distance < table.size(); ++distance)
        table[distance] = distance + 120u;

    for (unsigned code = 1; code <= kDistanceMap.size(); ++code) {
        const DistanceOffset offset = kDistanceMap[code - 1u];
        std::int64_t distance = static_cast<std::int64_t>(offset.x)
            + static_cast<std::int64_t>(offset.y) * imageWidth;
        if (distance < 1)
            distance = 1;
        if (distance >= static_cast<std::int64_t>(table.size()))
            continue;

        const unsigned scanDistance = static_cast<unsigned>(distance);
        const unsigned current = table[scanDistance];
        const unsigned candidateCost = distanceBitCost(code, distanceCodes);
        const unsigned currentCost = distanceBitCost(current, distanceCodes);
        if (candidateCost < currentCost
            || (candidateCost == currentCost && code < current))
            table[scanDistance] = code;
    }
    return table;
}

[[nodiscard]] Match findMatch(
    const std::vector<std::uint8_t>& rgba,
    std::size_t pixelCount,
    std::size_t position,
    const std::vector<std::int32_t>& heads,
    const std::vector<std::int32_t>& previous,
    std::size_t maxCandidates = 64) noexcept {
    constexpr std::size_t maxDistance = 65536;
    constexpr std::size_t maxLengthLimit = 4096;

    if (position + 2u >= pixelCount)
        return {};

    std::int32_t candidate = heads[hash3Pixels(rgba, position)];
    Match best;
    std::size_t searched = 0;
    const std::size_t maxLength = std::min(maxLengthLimit, pixelCount - position);

    while (candidate >= 0 && searched < maxCandidates) {
        const std::size_t candidatePosition = static_cast<std::size_t>(candidate);
        const std::size_t distance = position - candidatePosition;
        if (distance > maxDistance)
            break;

        if (pixelAt(rgba, candidatePosition) == pixelAt(rgba, position)
            && pixelAt(rgba, candidatePosition + 1u) == pixelAt(rgba, position + 1u)
            && pixelAt(rgba, candidatePosition + 2u) == pixelAt(rgba, position + 2u)) {
            std::size_t length = 3;
            while (length < maxLength
                && pixelAt(rgba, candidatePosition + length)
                    == pixelAt(rgba, position + length))
                ++length;
            if (length > best.length) {
                best = Match{length, distance};
                if (length == maxLength)
                    break;
            }
        }
        candidate = previous[candidatePosition];
        ++searched;
    }
    return best;
}

void insertPosition(
    const std::vector<std::uint8_t>& rgba,
    std::size_t pixelCount,
    std::size_t position,
    std::vector<std::int32_t>& heads,
    std::vector<std::int32_t>& previous) noexcept {
    if (position + 2u >= pixelCount)
        return;
    const std::uint16_t hash = hash3Pixels(rgba, position);
    previous[position] = heads[hash];
    heads[hash] = static_cast<std::int32_t>(position);
}

void rollbackPosition(
    const std::vector<std::uint8_t>& rgba,
    std::size_t pixelCount,
    std::size_t position,
    std::vector<std::int32_t>& heads,
    std::vector<std::int32_t>& previous) noexcept {
    if (position + 2u >= pixelCount)
        return;
    const std::uint16_t hash = hash3Pixels(rgba, position);
    if (heads[hash] == static_cast<std::int32_t>(position))
        heads[hash] = previous[position];
    previous[position] = -1;
}

[[nodiscard]] std::size_t literalPixelBits(
    const std::vector<std::uint8_t>& rgba,
    std::size_t position,
    const std::vector<PrefixCode>& greenCodes) noexcept {
    const std::uint8_t green = rgba[position * 4u + 1u];
    // R/B/Aはuniform 8-bit code，Gだけgreen/length treeを共有する。
    return static_cast<std::size_t>(greenCodes[green].bits) + 24u;
}

[[nodiscard]] std::size_t literalBits(
    const std::vector<std::uint8_t>& rgba,
    std::size_t position,
    std::size_t length,
    const std::vector<PrefixCode>& greenCodes) noexcept {
    std::size_t bits = 0;
    for (std::size_t i = 0; i < length; ++i)
        bits += literalPixelBits(rgba, position + i, greenCodes);
    return bits;
}

[[nodiscard]] std::size_t matchBits(
    const Match& match,
    const std::vector<unsigned>& mappedDistances,
    const std::vector<PrefixCode>& greenCodes,
    const std::vector<PrefixCode>& distanceCodes) noexcept {
    const PrefixInteger length = encodePrefixInteger(static_cast<unsigned>(match.length));
    const PrefixInteger distance = encodePrefixInteger(mappedDistances[match.distance]);
    return static_cast<std::size_t>(greenCodes[256u + length.code].bits)
        + length.extraBits
        + distanceCodes[distance.code].bits
        + distance.extraBits;
}

[[nodiscard]] std::size_t tokenSpan(const Match& match) noexcept {
    return match.length >= 3u ? match.length : 1u;
}

[[nodiscard]] std::size_t tokenBits(
    const std::vector<std::uint8_t>& rgba,
    std::size_t position,
    const Match& match,
    const std::vector<unsigned>& mappedDistances,
    const std::vector<PrefixCode>& greenCodes,
    const std::vector<PrefixCode>& distanceCodes) noexcept {
    if (match.length >= 3u)
        return matchBits(match, mappedDistances, greenCodes, distanceCodes);
    return literalPixelBits(rgba, position, greenCodes);
}

enum class ImageTokenKind : std::uint8_t {
    Literal,
    Match,
    ColorCache
};

struct ImageToken final {
    ImageTokenKind kind = ImageTokenKind::Literal;
    std::uint8_t red = 0u;
    std::uint8_t green = 0u;
    std::uint8_t blue = 0u;
    std::uint8_t alpha = 0u;
    std::uint16_t cacheIndex = 0u;
    std::uint16_t length = 0u;
    std::uint32_t distance = 0u; // scan-line pixel distance before VP8L 2D mapping
};

[[nodiscard]] std::vector<ImageToken> tokenizeImageData(
    const std::vector<std::uint8_t>& rgba,
    unsigned width,
    unsigned height,
    bool useColorCache,
    const std::vector<PrefixCode>& provisionalGreenCodes,
    const std::vector<PrefixCode>& provisionalDistanceCodes,
    const std::vector<unsigned>& mappedDistances) {
    const std::size_t pixelCount = static_cast<std::size_t>(width) * height;
    std::vector<ImageToken> tokens;
    tokens.reserve(std::min<std::size_t>(pixelCount, 65536u));

    std::vector<std::int32_t> heads(65536u, -1);
    std::vector<std::int32_t> previous(pixelCount, -1);
    std::optional<ColorCache> colorCache;
    if (useColorCache)
        colorCache.emplace(4u);

    const auto emitPixel = [&](std::size_t position) {
        const std::size_t i = position * 4u;
        if (colorCache) {
            const std::uint32_t argb = webpArgbAt(rgba, position);
            if (const auto cacheIndex = colorCache->find(argb)) {
                ImageToken token;
                token.kind = ImageTokenKind::ColorCache;
                token.cacheIndex = static_cast<std::uint16_t>(*cacheIndex);
                tokens.push_back(token);
                colorCache->insert(argb);
                return;
            }
        }

        ImageToken token;
        token.kind = ImageTokenKind::Literal;
        token.red = rgba[i + 0u];
        token.green = rgba[i + 1u];
        token.blue = rgba[i + 2u];
        token.alpha = rgba[i + 3u];
        tokens.push_back(token);
        if (colorCache)
            colorCache->insert(webpArgbAt(rgba, position));
    };

    std::size_t position = 0u;
    while (position < pixelCount) {
        const Match match = findMatch(rgba, pixelCount, position, heads, previous);

        bool wholeMatchInserted = false;
        if (match.length >= 3u) {
            // LZ77 parser自体はHuffman化前と同一に保つ。
            // provisional treeで従来の1-pixel lazy判定を行い，entropy treeだけを後段で最適化する。
            constexpr std::size_t maxLazyLength = 32u;
            constexpr std::size_t maxLazyCandidates = 16u;
            if (match.length <= maxLazyLength
                && position + match.length < pixelCount) {
                insertPosition(rgba, pixelCount, position, heads, previous);
                const Match lazyNext = findMatch(
                    rgba, pixelCount, position + 1u, heads, previous, maxLazyCandidates);

                for (std::size_t i = 1u; i < match.length; ++i)
                    insertPosition(rgba, pixelCount, position + i, heads, previous);
                wholeMatchInserted = true;
                const Match greedyNext = findMatch(
                    rgba, pixelCount, position + match.length,
                    heads, previous, maxLazyCandidates);

                const std::size_t lazySpan = 1u + tokenSpan(lazyNext);
                const std::size_t greedySpan = match.length + tokenSpan(greedyNext);
                const std::size_t lazyCost = literalPixelBits(
                        rgba, position, provisionalGreenCodes)
                    + tokenBits(
                        rgba, position + 1u, lazyNext,
                        mappedDistances, provisionalGreenCodes, provisionalDistanceCodes);
                const std::size_t greedyCost = matchBits(
                        match, mappedDistances,
                        provisionalGreenCodes, provisionalDistanceCodes)
                    + tokenBits(
                        rgba, position + match.length, greedyNext,
                        mappedDistances, provisionalGreenCodes, provisionalDistanceCodes);
                const std::ptrdiff_t lazySavings =
                    static_cast<std::ptrdiff_t>(literalBits(
                        rgba, position, lazySpan, provisionalGreenCodes))
                    - static_cast<std::ptrdiff_t>(lazyCost);
                const std::ptrdiff_t greedySavings =
                    static_cast<std::ptrdiff_t>(literalBits(
                        rgba, position, greedySpan, provisionalGreenCodes))
                    - static_cast<std::ptrdiff_t>(greedyCost);

                const bool chooseLazy =
                    (lazySpan >= greedySpan + 2u && lazySavings > greedySavings)
                    || (lazySpan == greedySpan && lazyCost < greedyCost);
                if (chooseLazy) {
                    for (std::size_t i = match.length; i-- > 1u;)
                        rollbackPosition(rgba, pixelCount, position + i, heads, previous);
                    emitPixel(position);
                    ++position;
                    continue;
                }
            }

            ImageToken token;
            token.kind = ImageTokenKind::Match;
            token.length = static_cast<std::uint16_t>(match.length);
            token.distance = static_cast<std::uint32_t>(match.distance);
            tokens.push_back(token);

            if (!wholeMatchInserted) {
                for (std::size_t i = 0; i < match.length; ++i)
                    insertPosition(rgba, pixelCount, position + i, heads, previous);
            }
            if (colorCache) {
                for (std::size_t i = 0; i < match.length; ++i)
                    colorCache->insert(webpArgbAt(rgba, position + i));
            }
            position += match.length;
            continue;
        }

        emitPixel(position);
        insertPosition(rgba, pixelCount, position, heads, previous);
        ++position;
    }
    return tokens;
}

struct ImageFrequencies final {
    std::vector<std::uint64_t> green;
    std::vector<std::uint64_t> red = std::vector<std::uint64_t>(256u, 0u);
    std::vector<std::uint64_t> blue = std::vector<std::uint64_t>(256u, 0u);
    std::vector<std::uint64_t> alpha = std::vector<std::uint64_t>(256u, 0u);
    std::vector<std::uint64_t> distance = std::vector<std::uint64_t>(40u, 0u);
};

[[nodiscard]] ImageFrequencies collectImageFrequencies(
    const std::vector<ImageToken>& tokens,
    unsigned colorCacheSize,
    const std::vector<unsigned>& mappedDistances) {
    ImageFrequencies frequencies;
    frequencies.green.assign(280u + colorCacheSize, 0u);

    for (const ImageToken& token : tokens) {
        switch (token.kind) {
        case ImageTokenKind::Literal:
            ++frequencies.green[token.green];
            ++frequencies.red[token.red];
            ++frequencies.blue[token.blue];
            ++frequencies.alpha[token.alpha];
            break;
        case ImageTokenKind::ColorCache:
            ++frequencies.green[280u + token.cacheIndex];
            break;
        case ImageTokenKind::Match: {
            const PrefixInteger length = encodePrefixInteger(token.length);
            ++frequencies.green[256u + length.code];
            const unsigned mappedDistance = mappedDistances[token.distance];
            const PrefixInteger distance = encodePrefixInteger(mappedDistance);
            ++frequencies.distance[distance.code];
            break;
        }
        }
    }
    return frequencies;
}

struct ImagePrefixTrees final {
    std::vector<std::uint8_t> greenLengths;
    std::vector<std::uint8_t> redLengths;
    std::vector<std::uint8_t> blueLengths;
    std::vector<std::uint8_t> alphaLengths;
    std::vector<std::uint8_t> distanceLengths;
    std::vector<PrefixCode> greenCodes;
    std::vector<PrefixCode> redCodes;
    std::vector<PrefixCode> blueCodes;
    std::vector<PrefixCode> alphaCodes;
    std::vector<PrefixCode> distanceCodes;
};

[[nodiscard]] std::optional<ImagePrefixTrees> makeImagePrefixTrees(
    const ImageFrequencies& frequencies) {
    ImagePrefixTrees trees;
    trees.greenLengths = makeHuffmanLengths(frequencies.green, 15u);
    trees.redLengths = makeHuffmanLengths(frequencies.red, 15u);
    trees.blueLengths = makeHuffmanLengths(frequencies.blue, 15u);
    trees.alphaLengths = makeHuffmanLengths(frequencies.alpha, 15u);
    trees.distanceLengths = makeHuffmanLengths(frequencies.distance, 15u);
    if (trees.greenLengths.empty() || trees.redLengths.empty()
        || trees.blueLengths.empty() || trees.alphaLengths.empty()
        || trees.distanceLengths.empty())
        return std::nullopt;

    trees.greenCodes = makeCanonicalCodes(trees.greenLengths);
    trees.redCodes = makeCanonicalCodes(trees.redLengths);
    trees.blueCodes = makeCanonicalCodes(trees.blueLengths);
    trees.alphaCodes = makeCanonicalCodes(trees.alphaLengths);
    trees.distanceCodes = makeCanonicalCodes(trees.distanceLengths);
    if (trees.greenCodes.size() != trees.greenLengths.size()
        || trees.redCodes.size() != trees.redLengths.size()
        || trees.blueCodes.size() != trees.blueLengths.size()
        || trees.alphaCodes.size() != trees.alphaLengths.size()
        || trees.distanceCodes.size() != trees.distanceLengths.size())
        return std::nullopt;
    return trees;
}

[[nodiscard]] bool writeEncodedImageData(
    BitWriter& writer,
    const std::vector<std::uint8_t>& rgba,
    unsigned width,
    unsigned height,
    bool spatiallyCoded,
    bool useColorCache) {
    const std::uint64_t pixelCount64 = static_cast<std::uint64_t>(width) * height;
    if (width == 0u || height == 0u
        || pixelCount64 > static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())
        || pixelCount64 > std::numeric_limits<std::size_t>::max() / 4u
        || rgba.size() != static_cast<std::size_t>(pixelCount64 * 4u))
        return false;

    constexpr unsigned colorCacheBits = 4u;
    const unsigned colorCacheSize = useColorCache ? (1u << colorCacheBits) : 0u;

    // LZ77/lazy parserは従来treeで固定し，今回の変更をentropy codingだけに限定する。
    const std::vector<std::uint8_t> provisionalGreenLengths =
        makeProvisionalGreenLengths(colorCacheSize);
    const auto provisionalGreenCodes = makeCanonicalCodes(provisionalGreenLengths);
    std::vector<std::uint8_t> provisionalDistanceLengths(
        40u, static_cast<std::uint8_t>(6));
    for (unsigned symbol = 0; symbol < 24u; ++symbol)
        provisionalDistanceLengths[symbol] = static_cast<std::uint8_t>(5);
    const auto provisionalDistanceCodes = makeCanonicalCodes(provisionalDistanceLengths);
    if (provisionalGreenCodes.size() != provisionalGreenLengths.size()
        || provisionalDistanceCodes.size() != provisionalDistanceLengths.size())
        return false;
    const auto mappedDistances = makeMappedDistanceTable(width, provisionalDistanceCodes);

    const std::vector<ImageToken> tokens = tokenizeImageData(
        rgba,
        width,
        height,
        useColorCache,
        provisionalGreenCodes,
        provisionalDistanceCodes,
        mappedDistances);
    const ImageFrequencies frequencies = collectImageFrequencies(
        tokens, colorCacheSize, mappedDistances);
    const auto trees = makeImagePrefixTrees(frequencies);
    if (!trees)
        return false;

    if (useColorCache) {
        writer.writeBits(1u, 1);
        writer.writeBits(colorCacheBits, 4);
    }
    else {
        writer.writeBits(0u, 1);
    }
    if (spatiallyCoded)
        writer.writeBits(0u, 1); // one meta prefix code group

    // VP8L prefix-code group order: green/length/cache, red, blue, alpha, distance.
    writePrefixCodeLengths(writer, trees->greenLengths);
    writePrefixCodeLengths(writer, trees->redLengths);
    writePrefixCodeLengths(writer, trees->blueLengths);
    writePrefixCodeLengths(writer, trees->alphaLengths);
    writePrefixCodeLengths(writer, trees->distanceLengths);

    for (const ImageToken& token : tokens) {
        switch (token.kind) {
        case ImageTokenKind::Literal:
            writePrefixSymbol(writer, trees->greenCodes, token.green);
            writePrefixSymbol(writer, trees->redCodes, token.red);
            writePrefixSymbol(writer, trees->blueCodes, token.blue);
            writePrefixSymbol(writer, trees->alphaCodes, token.alpha);
            break;

        case ImageTokenKind::ColorCache:
            writePrefixSymbol(writer, trees->greenCodes, 280u + token.cacheIndex);
            break;

        case ImageTokenKind::Match: {
            const PrefixInteger length = encodePrefixInteger(token.length);
            writePrefixSymbol(writer, trees->greenCodes, 256u + length.code);
            if (length.extraBits != 0u)
                writer.writeBits(length.extraValue, length.extraBits);

            const unsigned mappedDistance = mappedDistances[token.distance];
            const PrefixInteger distance = encodePrefixInteger(mappedDistance);
            writePrefixSymbol(writer, trees->distanceCodes, distance.code);
            if (distance.extraBits != 0u)
                writer.writeBits(distance.extraValue, distance.extraBits);
            break;
        }
        }
    }
    return true;
}

enum class WebpTransformPlan {
    None,
    Predictor,
    SubtractGreen,
    SubtractGreenPredictor
};

[[nodiscard]] RasterImage subtractGreen(const RasterImage& image) {
    RasterImage transformed = image;
    for (std::size_t i = 0; i < transformed.rgba.size(); i += 4u) {
        const std::uint8_t green = transformed.rgba[i + 1u];
        transformed.rgba[i + 0u] = static_cast<std::uint8_t>(
            static_cast<unsigned>(transformed.rgba[i + 0u])
            - static_cast<unsigned>(green));
        transformed.rgba[i + 2u] = static_cast<std::uint8_t>(
            static_cast<unsigned>(transformed.rgba[i + 2u])
            - static_cast<unsigned>(green));
    }
    return transformed;
}

[[nodiscard]] std::optional<std::string> encodeVp8lPayload(
    const RasterImage& image,
    bool alphaUsed,
    WebpTransformPlan plan) {
    BitWriter writer;
    writer.writeBits(image.widthPx - 1u, 14);
    writer.writeBits(image.heightPx - 1u, 14);
    writer.writeBits(alphaUsed ? 1u : 0u, 1);
    writer.writeBits(0u, 3); // version

    const bool useSubtractGreen = plan == WebpTransformPlan::SubtractGreen
        || plan == WebpTransformPlan::SubtractGreenPredictor;
    const bool usePredictor = plan == WebpTransformPlan::Predictor
        || plan == WebpTransformPlan::SubtractGreenPredictor;

    std::optional<RasterImage> subtractGreenImage;
    const RasterImage* transformedImage = &image;
    if (useSubtractGreen) {
        subtractGreenImage = subtractGreen(image);
        transformedImage = &*subtractGreenImage;

        // Subtract Green Transformには付随dataがない。forward transformでは
        // R'=R-G, B'=B-G (mod 256) とし，decoderがGを足して復元する。
        writer.writeBits(1u, 1); // transform present
        writer.writeBits(2u, 2); // subtract green transform
    }

    if (usePredictor) {
        // Subtract Greenと併用する場合は，subtract後のpixel列へPredictorを掛ける。
        // bitstreamにもSubtract Green -> Predictorの順で記録するため，decoderは
        // Predictor^-1 -> Add Greenの逆順で原画像へ戻す。
        constexpr unsigned predictorSizeBits = 4u;
        const PredictorTransform predictor = makePredictorTransform(
            *transformedImage, predictorSizeBits);

        writer.writeBits(1u, 1); // transform present
        writer.writeBits(0u, 2); // predictor transform
        writer.writeBits(predictorSizeBits - 2u, 3);
        // Predictor imageはentropy-coded subresolution imageなのでtransform terminatorも
        // meta-prefix bitも持たない。modeはgreen channel 0..13に格納する。
        if (!writeEncodedImageData(
                writer,
                predictor.modesRgba,
                predictor.transformWidth,
                predictor.transformHeight,
                false,
                false))
            return std::nullopt;

        writer.writeBits(0u, 1); // end of transforms
        if (!writeEncodedImageData(
                writer,
                predictor.residualRgba,
                image.widthPx,
                image.heightPx,
                true,
                true))
            return std::nullopt;
    }
    else {
        writer.writeBits(0u, 1); // end of transforms
        if (!writeEncodedImageData(
                writer,
                transformedImage->rgba,
                image.widthPx,
                image.heightPx,
                true,
                true))
            return std::nullopt;
    }

    std::vector<std::uint8_t> bits = writer.finish();
    std::string vp8l;
    vp8l.reserve(bits.size() + 1u);
    vp8l.push_back(static_cast<char>(0x2f));
    vp8l.append(reinterpret_cast<const char*>(bits.data()), bits.size());
    if (vp8l.size() > std::numeric_limits<std::uint32_t>::max())
        return std::nullopt;
    return vp8l;
}

[[nodiscard]] std::optional<std::string> encodeVp8l(const RasterImage& image) {
    if (image.widthPx == 0 || image.heightPx == 0
        || image.widthPx > 16384u || image.heightPx > 16384u)
        return std::nullopt;
    const std::uint64_t pixelCount64 =
        static_cast<std::uint64_t>(image.widthPx) * image.heightPx;
    if (pixelCount64 > std::numeric_limits<std::size_t>::max() / 4u
        || pixelCount64 > static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())
        || image.rgba.size() != static_cast<std::size_t>(pixelCount64 * 4u))
        return std::nullopt;

    bool alphaUsed = false;
    for (std::size_t i = 3; i < image.rgba.size(); i += 4) {
        if (image.rgba[i] != 255u) {
            alphaUsed = true;
            break;
        }
    }

    auto plain = encodeVp8lPayload(
        image, alphaUsed, WebpTransformPlan::None);
    auto predictor = encodeVp8lPayload(
        image, alphaUsed, WebpTransformPlan::Predictor);
    auto subtract = encodeVp8lPayload(
        image, alphaUsed, WebpTransformPlan::SubtractGreen);
    auto subtractPredictor = encodeVp8lPayload(
        image, alphaUsed, WebpTransformPlan::SubtractGreenPredictor);
    if (!plain || !predictor || !subtract || !subtractPredictor)
        return std::nullopt;

    // Transformは画像によってLZ77/color-cacheとの相性が逆転するため，
    // None / Predictor / Subtract Green / Subtract Green+Predictorを実際にencodeして
    // 最小payloadだけ採用する。同sizeでは単純な候補を優先する。
    std::string vp8l = std::move(*plain);
    const auto chooseSmaller = [&](std::optional<std::string>& candidate) {
        if (candidate->size() < vp8l.size())
            vp8l = std::move(*candidate);
    };
    chooseSmaller(predictor);
    chooseSmaller(subtract);
    chooseSmaller(subtractPredictor);

    std::string output;
    output.reserve(20u + vp8l.size() + (vp8l.size() & 1u));
    output.append("RIFF", 4);
    appendLittleEndian32(output, 0); // RIFF sizeは最後にpatchする。
    output.append("WEBP", 4);
    output.append("VP8L", 4);
    appendLittleEndian32(output, static_cast<std::uint32_t>(vp8l.size()));
    output.append(vp8l);
    if ((vp8l.size() & 1u) != 0u)
        output.push_back('\0'); // RIFF chunk padding。VP8L chunk sizeには含めない。

    if (output.size() < 8u || output.size() - 8u > std::numeric_limits<std::uint32_t>::max())
        return std::nullopt;
    const std::uint32_t riffSize = static_cast<std::uint32_t>(output.size() - 8u);
    output[4] = static_cast<char>(riffSize & 0xffu);
    output[5] = static_cast<char>((riffSize >> 8u) & 0xffu);
    output[6] = static_cast<char>((riffSize >> 16u) & 0xffu);
    output[7] = static_cast<char>((riffSize >> 24u) & 0xffu);
    return output;
}

} // namespace

WebpRenderResult renderWebp(const GraphicsScene& scene, const RasterRenderOptions& options) {
    const auto raster = renderRaster(scene, options);
    if (!raster) {
        switch (raster.status) {
        case RasterRenderStatus::InvalidScene:
            return WebpRenderResult{WebpRenderStatus::InvalidScene, std::nullopt};
        case RasterRenderStatus::InvalidOptions:
            return WebpRenderResult{WebpRenderStatus::InvalidOptions, std::nullopt};
        case RasterRenderStatus::UnsupportedText:
            return WebpRenderResult{WebpRenderStatus::UnsupportedText, std::nullopt};
        case RasterRenderStatus::ResourceLimit:
            return WebpRenderResult{WebpRenderStatus::ResourceLimit, std::nullopt};
        default:
            return WebpRenderResult{WebpRenderStatus::EncodeFailed, std::nullopt};
        }
    }

    try {
        auto encoded = encodeVp8l(*raster.image);
        if (!encoded)
            return WebpRenderResult{WebpRenderStatus::EncodeFailed, std::nullopt};
        return WebpRenderResult{WebpRenderStatus::Success, std::move(encoded)};
    }
    catch (const std::bad_alloc&) {
        return WebpRenderResult{WebpRenderStatus::ResourceLimit, std::nullopt};
    }
    catch (...) {
        return WebpRenderResult{WebpRenderStatus::EncodeFailed, std::nullopt};
    }
}

} // namespace mmcal::graphics
