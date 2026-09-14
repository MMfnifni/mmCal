#include "png_backend.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <queue>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mmcal::graphics {
namespace {

[[nodiscard]] std::uint32_t crc32(std::string_view data) noexcept {
    std::uint32_t crc = 0xffffffffu;
    for (const unsigned char byte : data) {
        crc ^= byte;
        for (int bit = 0; bit < 8; ++bit)
            crc = (crc >> 1u) ^ (0xedb88320u & (0u - (crc & 1u)));
    }
    return ~crc;
}

[[nodiscard]] std::uint32_t adler32(const std::vector<std::uint8_t>& data) noexcept {
    constexpr std::uint32_t modulus = 65521u;
    std::uint32_t a = 1u;
    std::uint32_t b = 0u;
    std::size_t offset = 0;
    while (offset < data.size()) {
        const std::size_t end = std::min(data.size(), offset + 5552u);
        for (; offset < end; ++offset) {
            a += data[offset];
            b += a;
        }
        a %= modulus;
        b %= modulus;
    }
    return (b << 16u) | a;
}

void appendBigEndian32(std::string& output, std::uint32_t value) {
    output.push_back(static_cast<char>((value >> 24u) & 0xffu));
    output.push_back(static_cast<char>((value >> 16u) & 0xffu));
    output.push_back(static_cast<char>((value >> 8u) & 0xffu));
    output.push_back(static_cast<char>(value & 0xffu));
}

void appendChunk(std::string& output, std::string_view type, std::string_view data) {
    appendBigEndian32(output, static_cast<std::uint32_t>(data.size()));
    const std::size_t crcStart = output.size();
    output.append(type);
    output.append(data);
    appendBigEndian32(output, crc32(std::string_view{output}.substr(crcStart, 4u + data.size())));
}

[[nodiscard]] std::uint8_t paethPredictor(
    std::uint8_t left,
    std::uint8_t up,
    std::uint8_t upperLeft) noexcept {
    const int p = static_cast<int>(left) + static_cast<int>(up) - static_cast<int>(upperLeft);
    const int pa = std::abs(p - static_cast<int>(left));
    const int pb = std::abs(p - static_cast<int>(up));
    const int pc = std::abs(p - static_cast<int>(upperLeft));
    if (pa <= pb && pa <= pc)
        return left;
    if (pb <= pc)
        return up;
    return upperLeft;
}

[[nodiscard]] std::uint64_t filterScore(const std::vector<std::uint8_t>& row) noexcept {
    std::uint64_t score = 0;
    for (const std::uint8_t value : row)
        score += std::min<std::uint32_t>(value, 256u - value);
    return score;
}

[[nodiscard]] std::vector<std::uint8_t> filterScanlines(const RasterImage& image) {
    constexpr std::size_t bytesPerPixel = 4;
    const std::size_t rowBytes = static_cast<std::size_t>(image.widthPx) * bytesPerPixel;
    std::vector<std::uint8_t> filtered;
    filtered.reserve((rowBytes + 1u) * image.heightPx);
    std::array<std::vector<std::uint8_t>, 5> candidates;
    for (auto& candidate : candidates)
        candidate.resize(rowBytes);

    for (std::uint32_t y = 0; y < image.heightPx; ++y) {
        const std::uint8_t* current = image.rgba.data() + static_cast<std::size_t>(y) * rowBytes;
        const std::uint8_t* previous = y == 0
            ? nullptr
            : image.rgba.data() + static_cast<std::size_t>(y - 1) * rowBytes;
        for (std::size_t x = 0; x < rowBytes; ++x) {
            const std::uint8_t left = x >= bytesPerPixel ? current[x - bytesPerPixel] : 0;
            const std::uint8_t up = previous ? previous[x] : 0;
            const std::uint8_t upperLeft = previous && x >= bytesPerPixel
                ? previous[x - bytesPerPixel]
                : 0;
            candidates[0][x] = current[x];
            candidates[1][x] = static_cast<std::uint8_t>(current[x] - left);
            candidates[2][x] = static_cast<std::uint8_t>(current[x] - up);
            candidates[3][x] = static_cast<std::uint8_t>(
                current[x] - static_cast<std::uint8_t>(
                    (static_cast<unsigned>(left) + static_cast<unsigned>(up)) / 2u));
            candidates[4][x] = static_cast<std::uint8_t>(
                current[x] - paethPredictor(left, up, upperLeft));
        }

        std::size_t best = 0;
        std::uint64_t bestScore = filterScore(candidates[0]);
        for (std::size_t filter = 1; filter < candidates.size(); ++filter) {
            const std::uint64_t score = filterScore(candidates[filter]);
            if (score < bestScore) {
                best = filter;
                bestScore = score;
            }
        }
        filtered.push_back(static_cast<std::uint8_t>(best));
        filtered.insert(filtered.end(), candidates[best].begin(), candidates[best].end());
    }
    return filtered;
}

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

struct HuffmanCode final {
    std::uint16_t code = 0;
    std::uint8_t bits = 0;
};

// 頻度からDEFLATE用canonical Huffman codeのbit-lengthを作る。
// data treeは最大15 bit，code-length treeは最大7 bit。
// 通常木が上限を越える極端な分布では，完全木のbalanced fallbackへ落とす。
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

    // DEFLATEでもtreeは最低1 symbol必要。empty distance treeはsymbol 0を使う。
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

[[nodiscard]] std::vector<HuffmanCode> makeCanonicalCodes(
    const std::vector<std::uint8_t>& lengths,
    unsigned maxBits) {
    std::vector<unsigned> count(maxBits + 1u, 0u);
    for (const std::uint8_t length : lengths) {
        if (length > maxBits)
            return {};
        if (length != 0u)
            ++count[length];
    }

    std::vector<unsigned> next(maxBits + 1u, 0u);
    unsigned code = 0u;
    for (unsigned bits = 1u; bits <= maxBits; ++bits) {
        code = (code + count[bits - 1u]) << 1u;
        next[bits] = code;
    }

    std::vector<HuffmanCode> table(lengths.size());
    for (std::size_t symbol = 0; symbol < lengths.size(); ++symbol) {
        const unsigned bits = lengths[symbol];
        if (bits == 0u)
            continue;
        table[symbol] = HuffmanCode{
            static_cast<std::uint16_t>(next[bits]++),
            static_cast<std::uint8_t>(bits)};
    }
    return table;
}

void writeHuffmanSymbol(
    BitWriter& writer,
    const std::vector<HuffmanCode>& table,
    unsigned symbol) {
    const HuffmanCode code = table[symbol];
    writer.writeBits(reverseBits(code.code, code.bits), code.bits);
}

void writeFixedLiteralLength(BitWriter& writer, unsigned symbol) {
    std::uint32_t code = 0;
    unsigned bits = 0;
    if (symbol <= 143u) {
        code = 0x30u + symbol;
        bits = 8;
    }
    else if (symbol <= 255u) {
        code = 0x190u + (symbol - 144u);
        bits = 9;
    }
    else if (symbol <= 279u) {
        code = symbol - 256u;
        bits = 7;
    }
    else {
        code = 0xc0u + (symbol - 280u);
        bits = 8;
    }
    writer.writeBits(reverseBits(code, bits), bits);
}

void writeFixedDistance(BitWriter& writer, unsigned symbol) {
    writer.writeBits(reverseBits(symbol, 5), 5);
}

struct LengthCode final {
    unsigned symbol = 257;
    unsigned extraBits = 0;
    unsigned extraValue = 0;
};

[[nodiscard]] LengthCode encodeLength(std::size_t length) {
    // RFC 1951ではcode 284は227..257，258は専用のcode 285で表す。
    // 284 + extra=31は復号器によっては258として受理されるが，仕様外で5 bit余計に掛かる。
    if (length >= 258u)
        return LengthCode{285u, 0u, 0u};

    static constexpr std::array<unsigned, 29> bases{
        3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,
        131,163,195,227,258};
    static constexpr std::array<unsigned, 29> extras{
        0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0};
    for (std::size_t i = 0; i < bases.size(); ++i) {
        const unsigned span = extras[i] == 0 ? 1u : (1u << extras[i]);
        if (length >= bases[i] && length < static_cast<std::size_t>(bases[i] + span))
            return LengthCode{
                static_cast<unsigned>(257u + i), extras[i],
                static_cast<unsigned>(length - bases[i])};
    }
    return LengthCode{285, 0, 0};
}

struct DistanceCode final {
    unsigned symbol = 0;
    unsigned extraBits = 0;
    unsigned extraValue = 0;
};

[[nodiscard]] DistanceCode encodeDistance(std::size_t distance) {
    static constexpr std::array<unsigned, 30> bases{
        1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,
        1537,2049,3073,4097,6145,8193,12289,16385,24577};
    static constexpr std::array<unsigned, 30> extras{
        0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13};
    for (std::size_t i = 0; i < bases.size(); ++i) {
        const unsigned span = extras[i] == 0 ? 1u : (1u << extras[i]);
        if (distance >= bases[i] && distance < static_cast<std::size_t>(bases[i] + span))
            return DistanceCode{
                static_cast<unsigned>(i), extras[i],
                static_cast<unsigned>(distance - bases[i])};
    }
    return DistanceCode{29, 13, static_cast<unsigned>(distance - 24577)};
}

[[nodiscard]] std::uint16_t matchHash(
    const std::vector<std::uint8_t>& input,
    std::size_t position,
    bool extended) noexcept {
    std::uint32_t value = static_cast<std::uint32_t>(input[position]) * 251u;
    value = (value ^ input[position + 1u]) * 251u;
    value ^= input[position + 2u];
    // PNG filter後は0が多く3-byte hashが衝突しやすい。一方，短いmatchでは
    // 3-byte chainが有利なので，3-byte/8-byteの2 parserを独立に試して最終byte数で選ぶ。
    if (extended && position + 7u < input.size()) {
        for (std::size_t i = 3u; i < 8u; ++i)
            value = (value ^ input[position + i]) * 251u;
    }
    value ^= value >> 16u;
    return static_cast<std::uint16_t>(value & 0xffffu);
}

struct Match final {
    std::size_t length = 0;
    std::size_t distance = 0;
};

[[nodiscard]] unsigned fixedLiteralLengthBits(unsigned symbol) noexcept {
    if (symbol <= 143u)
        return 8u;
    if (symbol <= 255u)
        return 9u;
    if (symbol <= 279u)
        return 7u;
    return 8u;
}

struct DeflateCostModel final {
    std::vector<std::uint8_t> literalLengthBits;
    std::vector<std::uint8_t> distanceBits;
};

[[nodiscard]] DeflateCostModel fixedDeflateCostModel() {
    DeflateCostModel model;
    model.literalLengthBits.resize(286u, static_cast<std::uint8_t>(0));
    for (unsigned symbol = 0; symbol < model.literalLengthBits.size(); ++symbol)
        model.literalLengthBits[symbol] = static_cast<std::uint8_t>(fixedLiteralLengthBits(symbol));
    model.distanceBits.resize(30u, static_cast<std::uint8_t>(5));
    return model;
}

[[nodiscard]] unsigned literalCodeBits(
    const DeflateCostModel& model,
    unsigned symbol) noexcept {
    if (symbol < model.literalLengthBits.size() && model.literalLengthBits[symbol] != 0u)
        return model.literalLengthBits[symbol];
    return fixedLiteralLengthBits(symbol);
}

[[nodiscard]] unsigned distanceCodeBits(
    const DeflateCostModel& model,
    unsigned symbol) noexcept {
    if (symbol < model.distanceBits.size() && model.distanceBits[symbol] != 0u)
        return model.distanceBits[symbol];
    return 5u;
}

[[nodiscard]] std::size_t literalBits(
    const std::vector<std::uint8_t>& input,
    std::size_t position,
    std::size_t length,
    const DeflateCostModel& model) noexcept {
    std::size_t bits = 0;
    const std::size_t end = position + std::min(length, input.size() - position);
    for (std::size_t i = position; i < end; ++i)
        bits += literalCodeBits(model, input[i]);
    return bits;
}

[[nodiscard]] std::size_t matchBits(
    const Match& match,
    const DeflateCostModel& model) noexcept {
    const LengthCode lengthCode = encodeLength(match.length);
    const DistanceCode distanceCode = encodeDistance(match.distance);
    return literalCodeBits(model, lengthCode.symbol)
        + lengthCode.extraBits
        + distanceCodeBits(model, distanceCode.symbol)
        + distanceCode.extraBits;
}

[[nodiscard]] std::ptrdiff_t matchSavingsBits(
    const std::vector<std::uint8_t>& input,
    std::size_t position,
    const Match& match,
    const DeflateCostModel& model) noexcept {
    if (match.length < 3)
        return std::numeric_limits<std::ptrdiff_t>::min();
    return static_cast<std::ptrdiff_t>(literalBits(input, position, match.length, model))
        - static_cast<std::ptrdiff_t>(matchBits(match, model));
}

[[nodiscard]] Match findMatch(
    const std::vector<std::uint8_t>& input,
    std::size_t position,
    const std::array<std::int32_t, 65536>& heads,
    const std::vector<std::int32_t>& previous,
    const DeflateCostModel& model,
    bool extendedHash,
    std::size_t maxCandidates = 128u) noexcept {
    constexpr std::size_t windowSize = 32768u;
    constexpr std::size_t niceLength = 128u;
    if (position + 2u >= input.size())
        return {};

    const std::uint16_t hash = matchHash(input, position, extendedHash);
    std::int32_t candidate = heads[hash];
    Match best;
    std::size_t bestCost = std::numeric_limits<std::size_t>::max();
    std::size_t searched = 0u;
    const std::size_t maxLength = std::min<std::size_t>(258u, input.size() - position);

    while (candidate >= 0 && searched < maxCandidates) {
        const std::size_t candidatePosition = static_cast<std::size_t>(candidate);
        const std::size_t distance = position - candidatePosition;
        if (distance > windowSize)
            break;

        // 同一hashでも末尾byteが異なる候補は早く捨てる。best未確定時は3 byteだけ確認する。
        if (input[candidatePosition] == input[position]
            && input[candidatePosition + 1u] == input[position + 1u]
            && input[candidatePosition + 2u] == input[position + 2u]
            && (best.length < 3u || best.length >= maxLength
                || input[candidatePosition + best.length] == input[position + best.length])) {
            std::size_t length = 3u;
            while (length < maxLength
                && input[candidatePosition + length] == input[position + length])
                ++length;

            if (length >= 3u) {
                const Match candidateMatch{length, distance};
                const std::size_t cost = matchBits(candidateMatch, model);
                if (length > best.length || (length == best.length && cost < bestCost)) {
                    best = candidateMatch;
                    bestCost = cost;
                    if (length == maxLength || length >= niceLength)
                        break;
                }
            }
        }

        candidate = previous[candidatePosition];
        ++searched;

        // 既に十分長いmatchが得られたらchainを深追いしない。
        if (best.length >= 64u && searched >= 32u)
            break;
        if (best.length >= 32u && searched >= 64u)
            break;
    }
    return best;
}

void insertPosition(
    const std::vector<std::uint8_t>& input,
    std::size_t position,
    std::array<std::int32_t, 65536>& heads,
    std::vector<std::int32_t>& previous,
    bool extendedHash) noexcept {
    if (position + 2u >= input.size())
        return;
    const std::uint16_t hash = matchHash(input, position, extendedHash);
    previous[position] = heads[hash];
    heads[hash] = static_cast<std::int32_t>(position);
}

void rollbackPosition(
    const std::vector<std::uint8_t>& input,
    std::size_t position,
    std::array<std::int32_t, 65536>& heads,
    std::vector<std::int32_t>& previous,
    bool extendedHash) noexcept {
    if (position + 2u >= input.size())
        return;
    const std::uint16_t hash = matchHash(input, position, extendedHash);
    if (heads[hash] == static_cast<std::int32_t>(position))
        heads[hash] = previous[position];
    previous[position] = -1;
}

struct DeflateToken final {
    bool isMatch = false;
    std::uint8_t literal = 0u;
    Match match{};
};

[[nodiscard]] std::optional<DeflateCostModel> dynamicDeflateCostModel(
    const std::vector<DeflateToken>& tokens) {
    std::vector<std::uint64_t> literalLengthFreq(286u, 0u);
    std::vector<std::uint64_t> distanceFreq(30u, 0u);
    for (const DeflateToken& token : tokens) {
        if (!token.isMatch) {
            ++literalLengthFreq[token.literal];
            continue;
        }
        const LengthCode lengthCode = encodeLength(token.match.length);
        const DistanceCode distanceCode = encodeDistance(token.match.distance);
        ++literalLengthFreq[lengthCode.symbol];
        ++distanceFreq[distanceCode.symbol];
    }
    ++literalLengthFreq[256u];

    DeflateCostModel model;
    model.literalLengthBits = makeHuffmanLengths(literalLengthFreq, 15u);
    model.distanceBits = makeHuffmanLengths(distanceFreq, 15u);
    if (model.literalLengthBits.empty() || model.distanceBits.empty())
        return std::nullopt;
    return model;
}

struct EstimatedToken final {
    std::size_t length = 1u;
    std::size_t bits = 0u;
};

[[nodiscard]] EstimatedToken estimateNextToken(
    const std::vector<std::uint8_t>& input,
    std::size_t position,
    const std::array<std::int32_t, 65536>& heads,
    const std::vector<std::int32_t>& previous,
    const DeflateCostModel& model,
    bool extendedHash,
    std::size_t maxCandidates) noexcept {
    if (position >= input.size())
        return EstimatedToken{0u, 0u};
    Match match = findMatch(input, position, heads, previous, model, extendedHash, maxCandidates);
    if (match.length >= 3u && matchSavingsBits(input, position, match, model) > 0)
        return EstimatedToken{match.length, matchBits(match, model)};
    return EstimatedToken{1u, literalCodeBits(model, input[position])};
}

[[nodiscard]] std::vector<DeflateToken> tokenizeDeflate(
    const std::vector<std::uint8_t>& input,
    const DeflateCostModel& model,
    bool extendedHash) {
    std::array<std::int32_t, 65536> heads{};
    heads.fill(-1);
    std::vector<std::int32_t> previous(input.size(), -1);
    std::vector<DeflateToken> tokens;
    tokens.reserve(input.size() / 2u + 16u);

    std::size_t position = 0u;
    while (position < input.size()) {
        Match match = findMatch(input, position, heads, previous, model, extendedHash);
        if (match.length >= 3u && matchSavingsBits(input, position, match, model) <= 0)
            match = {};

        bool currentInserted = false;
        if (match.length >= 3u) {
            // 1-byte lazyを2-token局所評価へ拡張する。長いmatchはほぼ常に採用し，
            // 短中距離だけ current-match+次token と literal+次match+次token を比較する。
            constexpr std::size_t maxLazyLength = 96u;
            constexpr std::size_t lazyCandidates = 64u;
            constexpr std::size_t previewCandidates = 24u;
            if (match.length <= maxLazyLength && position + 1u < input.size()) {
                // greedy pathをpreviewするため，current match区間を一時的にhistoryへ入れる。
                for (std::size_t i = 0u; i < match.length; ++i)
                    insertPosition(input, position + i, heads, previous, extendedHash);
                const EstimatedToken greedyNext = estimateNextToken(
                    input, position + match.length, heads, previous, model, extendedHash, previewCandidates);
                const std::size_t greedyReach = match.length + greedyNext.length;
                const std::size_t greedyCost = matchBits(match, model) + greedyNext.bits;
                const std::size_t greedyLiteralCost = literalBits(
                    input, position, greedyReach, model);
                const std::ptrdiff_t greedySavings =
                    static_cast<std::ptrdiff_t>(greedyLiteralCost)
                    - static_cast<std::ptrdiff_t>(greedyCost);
                for (std::size_t i = match.length; i-- > 0u;)
                    rollbackPosition(input, position + i, heads, previous, extendedHash);

                // lazy path: current byteをliteralにして1 byte先のmatchを評価する。
                insertPosition(input, position, heads, previous, extendedHash);
                currentInserted = true;
                Match next = findMatch(
                    input, position + 1u, heads, previous, model, extendedHash, lazyCandidates);
                if (next.length >= 3u
                    && matchSavingsBits(input, position + 1u, next, model) > 0) {
                    for (std::size_t i = 0u; i < next.length; ++i)
                        insertPosition(input, position + 1u + i, heads, previous, extendedHash);
                    const EstimatedToken lazyNext = estimateNextToken(
                        input, position + 1u + next.length,
                        heads, previous, model, extendedHash, previewCandidates);
                    const std::size_t lazyReach = 1u + next.length + lazyNext.length;
                    const std::size_t lazyCost = literalCodeBits(model, input[position])
                        + matchBits(next, model) + lazyNext.bits;
                    const std::size_t lazyLiteralCost = literalBits(
                        input, position, lazyReach, model);
                    const std::ptrdiff_t lazySavings =
                        static_cast<std::ptrdiff_t>(lazyLiteralCost)
                        - static_cast<std::ptrdiff_t>(lazyCost);
                    for (std::size_t i = next.length; i-- > 0u;)
                        rollbackPosition(input, position + 1u + i, heads, previous, extendedHash);

                    if (lazyReach >= greedyReach && lazySavings > greedySavings) {
                        tokens.push_back(DeflateToken{false, input[position], {}});
                        ++position;
                        continue;
                    }
                }
            }

            tokens.push_back(DeflateToken{true, 0u, match});
            const std::size_t firstInsert = currentInserted ? 1u : 0u;
            for (std::size_t i = firstInsert; i < match.length; ++i)
                insertPosition(input, position + i, heads, previous, extendedHash);
            position += match.length;
        }
        else {
            tokens.push_back(DeflateToken{false, input[position], {}});
            insertPosition(input, position, heads, previous, extendedHash);
            ++position;
        }
    }
    return tokens;
}

void writeFixedTokens(BitWriter& writer, const std::vector<DeflateToken>& tokens) {
    for (const DeflateToken& token : tokens) {
        if (!token.isMatch) {
            writeFixedLiteralLength(writer, token.literal);
            continue;
        }
        const LengthCode lengthCode = encodeLength(token.match.length);
        const DistanceCode distanceCode = encodeDistance(token.match.distance);
        writeFixedLiteralLength(writer, lengthCode.symbol);
        if (lengthCode.extraBits != 0u)
            writer.writeBits(lengthCode.extraValue, lengthCode.extraBits);
        writeFixedDistance(writer, distanceCode.symbol);
        if (distanceCode.extraBits != 0u)
            writer.writeBits(distanceCode.extraValue, distanceCode.extraBits);
    }
    writeFixedLiteralLength(writer, 256u);
}

[[nodiscard]] std::vector<std::uint8_t> deflateFixed(
    const std::vector<DeflateToken>& tokens) {
    BitWriter writer;
    // BFINAL=1, BTYPE=01 (fixed Huffman). Bits are emitted least-significant first.
    writer.writeBits(1u, 1);
    writer.writeBits(1u, 2);
    writeFixedTokens(writer, tokens);
    return writer.finish();
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

[[nodiscard]] std::size_t lastNonzeroLength(
    const std::vector<std::uint8_t>& lengths,
    std::size_t minimumCount) noexcept {
    std::size_t count = std::min(minimumCount, lengths.size());
    for (std::size_t i = lengths.size(); i > count; --i) {
        if (lengths[i - 1u] != 0u)
            return i;
    }
    return count;
}

[[nodiscard]] std::vector<std::uint8_t> deflateDynamic(
    const std::vector<DeflateToken>& tokens) {
    std::vector<std::uint64_t> literalLengthFreq(286u, 0u);
    std::vector<std::uint64_t> distanceFreq(30u, 0u);
    for (const DeflateToken& token : tokens) {
        if (!token.isMatch) {
            ++literalLengthFreq[token.literal];
            continue;
        }
        const LengthCode lengthCode = encodeLength(token.match.length);
        const DistanceCode distanceCode = encodeDistance(token.match.distance);
        ++literalLengthFreq[lengthCode.symbol];
        ++distanceFreq[distanceCode.symbol];
    }
    ++literalLengthFreq[256u]; // end-of-block

    const std::vector<std::uint8_t> literalLengthLengths =
        makeHuffmanLengths(literalLengthFreq, 15u);
    const std::vector<std::uint8_t> distanceLengths =
        makeHuffmanLengths(distanceFreq, 15u);
    if (literalLengthLengths.empty() || distanceLengths.empty())
        return {};

    const std::size_t literalLengthCount =
        lastNonzeroLength(literalLengthLengths, 257u);
    const std::size_t distanceCount = lastNonzeroLength(distanceLengths, 1u);
    if (literalLengthCount > 286u || distanceCount > 30u)
        return {};

    std::vector<std::uint8_t> combinedLengths;
    combinedLengths.reserve(literalLengthCount + distanceCount);
    combinedLengths.insert(combinedLengths.end(),
        literalLengthLengths.begin(), literalLengthLengths.begin() + literalLengthCount);
    combinedLengths.insert(combinedLengths.end(),
        distanceLengths.begin(), distanceLengths.begin() + distanceCount);

    const std::vector<CodeLengthToken> codeLengthTokens =
        makeCodeLengthTokens(combinedLengths);
    std::vector<std::uint64_t> codeLengthFreq(19u, 0u);
    for (const CodeLengthToken& token : codeLengthTokens)
        ++codeLengthFreq[token.symbol];
    const std::vector<std::uint8_t> codeLengthLengths =
        makeHuffmanLengths(codeLengthFreq, 7u);
    if (codeLengthLengths.empty())
        return {};

    const auto literalLengthCodes = makeCanonicalCodes(literalLengthLengths, 15u);
    const auto distanceCodes = makeCanonicalCodes(distanceLengths, 15u);
    const auto codeLengthCodes = makeCanonicalCodes(codeLengthLengths, 7u);
    if (literalLengthCodes.empty() || distanceCodes.empty() || codeLengthCodes.empty())
        return {};

    static constexpr std::array<unsigned, 19> codeLengthOrder{
        16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15};
    unsigned codeLengthCount = 4u;
    for (unsigned i = 0; i < codeLengthOrder.size(); ++i) {
        if (codeLengthLengths[codeLengthOrder[i]] != 0u)
            codeLengthCount = std::max(codeLengthCount, i + 1u);
    }

    BitWriter writer;
    // BFINAL=1, BTYPE=10 (dynamic Huffman).
    writer.writeBits(1u, 1);
    writer.writeBits(2u, 2);
    writer.writeBits(static_cast<unsigned>(literalLengthCount - 257u), 5);
    writer.writeBits(static_cast<unsigned>(distanceCount - 1u), 5);
    writer.writeBits(codeLengthCount - 4u, 4);
    for (unsigned i = 0; i < codeLengthCount; ++i)
        writer.writeBits(codeLengthLengths[codeLengthOrder[i]], 3);

    for (const CodeLengthToken& token : codeLengthTokens) {
        writeHuffmanSymbol(writer, codeLengthCodes, token.symbol);
        if (token.extraBits != 0u)
            writer.writeBits(token.extraValue, token.extraBits);
    }

    for (const DeflateToken& token : tokens) {
        if (!token.isMatch) {
            writeHuffmanSymbol(writer, literalLengthCodes, token.literal);
            continue;
        }
        const LengthCode lengthCode = encodeLength(token.match.length);
        const DistanceCode distanceCode = encodeDistance(token.match.distance);
        writeHuffmanSymbol(writer, literalLengthCodes, lengthCode.symbol);
        if (lengthCode.extraBits != 0u)
            writer.writeBits(lengthCode.extraValue, lengthCode.extraBits);
        writeHuffmanSymbol(writer, distanceCodes, distanceCode.symbol);
        if (distanceCode.extraBits != 0u)
            writer.writeBits(distanceCode.extraValue, distanceCode.extraBits);
    }
    writeHuffmanSymbol(writer, literalLengthCodes, 256u);
    return writer.finish();
}

[[nodiscard]] std::string zlibDeflate(const std::vector<std::uint8_t>& input) {
    const DeflateCostModel fixedModel = fixedDeflateCostModel();
    std::vector<std::uint8_t> best;

    const auto consider = [&](std::vector<std::uint8_t> candidate) {
        if (!candidate.empty() && (best.empty() || candidate.size() < best.size()))
            best = std::move(candidate);
    };
    const auto runParser = [&](bool extendedHash) {
        const std::vector<DeflateToken> firstTokens =
            tokenizeDeflate(input, fixedModel, extendedHash);
        consider(deflateFixed(firstTokens));
        consider(deflateDynamic(firstTokens));

        // 実frequency treeを1回だけcost modelへ戻して再parseする。
        if (const auto dynamicModel = dynamicDeflateCostModel(firstTokens)) {
            const std::vector<DeflateToken> refinedTokens =
                tokenizeDeflate(input, *dynamicModel, extendedHash);
            consider(deflateDynamic(refinedTokens));
            consider(deflateFixed(refinedTokens));
        }
    };

    // 3-byte hashは短いmatchを拾いやすい基準parser。
    runParser(false);

    // 8-byte hashはPNG filter後の大量の3-byte衝突を避けられるが，parseをもう一巡要する。
    // 小中規模入力，または基準parserで圧縮率が低い入力だけ試し，bounded workを保つ。
    constexpr std::size_t extendedHashInputLimit = 1u << 20u;
    const bool weakCompression = !input.empty() && best.size() > input.size() / 20u;
    if (input.size() <= extendedHashInputLimit || weakCompression)
        runParser(true);

    const std::vector<std::uint8_t>& deflate = best;

    std::string output;
    output.reserve(deflate.size() + 6u);
    // CM=8, CINFO=7 (32 KiB window), FLEVEL=2. 0x789c is divisible by 31.
    output.push_back(static_cast<char>(0x78));
    output.push_back(static_cast<char>(0x9c));
    output.append(
        reinterpret_cast<const char*>(deflate.data()),
        static_cast<std::streamsize>(deflate.size()));
    appendBigEndian32(output, adler32(input));
    return output;
}

[[nodiscard]] std::string encodePng(const RasterImage& image) {
    const std::vector<std::uint8_t> filtered = filterScanlines(image);
    const std::string compressed = zlibDeflate(filtered);

    std::string png;
    png.reserve(compressed.size() + 128u);
    png.append("\x89PNG\r\n\x1a\n", 8);

    std::string ihdr;
    ihdr.reserve(13);
    appendBigEndian32(ihdr, image.widthPx);
    appendBigEndian32(ihdr, image.heightPx);
    ihdr.push_back(static_cast<char>(8));  // bit depth
    ihdr.push_back(static_cast<char>(6));  // RGBA
    ihdr.push_back(static_cast<char>(0));  // compression
    ihdr.push_back(static_cast<char>(0));  // filter
    ihdr.push_back(static_cast<char>(0));  // interlace
    appendChunk(png, "IHDR", ihdr);

    const double pixelsPerMeterX = image.dpiX / 0.0254;
    const double pixelsPerMeterY = image.dpiY / 0.0254;
    std::string phys;
    phys.reserve(9);
    appendBigEndian32(phys, static_cast<std::uint32_t>(std::llround(pixelsPerMeterX)));
    appendBigEndian32(phys, static_cast<std::uint32_t>(std::llround(pixelsPerMeterY)));
    phys.push_back(static_cast<char>(1));
    appendChunk(png, "pHYs", phys);

    appendChunk(png, "IDAT", compressed);
    appendChunk(png, "IEND", {});
    return png;
}

} // namespace

PngRenderResult renderPng(const GraphicsScene& scene, const RasterRenderOptions& options) {
    const auto raster = renderRaster(scene, options);
    if (!raster || !raster.image) {
        PngRenderStatus status = PngRenderStatus::EncodeFailed;
        switch (raster.status) {
        case RasterRenderStatus::InvalidScene: status = PngRenderStatus::InvalidScene; break;
        case RasterRenderStatus::InvalidOptions: status = PngRenderStatus::InvalidOptions; break;
        case RasterRenderStatus::UnsupportedText: status = PngRenderStatus::UnsupportedText; break;
        case RasterRenderStatus::ResourceLimit: status = PngRenderStatus::ResourceLimit; break;
        case RasterRenderStatus::RenderFailed: status = PngRenderStatus::EncodeFailed; break;
        case RasterRenderStatus::Success: status = PngRenderStatus::EncodeFailed; break;
        }
        return PngRenderResult{status, std::nullopt};
    }
    try {
        return PngRenderResult{PngRenderStatus::Success, encodePng(*raster.image)};
    }
    catch (...) {
        return PngRenderResult{PngRenderStatus::EncodeFailed, std::nullopt};
    }
}

} // namespace mmcal::graphics
