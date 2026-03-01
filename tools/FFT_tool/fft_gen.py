import random

def write_fft_random_matrix(filename, rows, cols, seed=1234):
    random.seed(seed)

    with open(filename, "w", encoding="utf-8") as f:
        f.write("{")

        for i in range(rows):
            f.write("{")
            for j in range(cols):
                val = random.uniform(-1.0, 1.0)
                f.write(f"{val:.10f}")
                if j != cols - 1:
                    f.write(",")
            f.write("}")
            if i != rows - 1:
                f.write(",")

        f.write("}")

    print(f"Written: {filename}")


write_fft_random_matrix("fft_random_1536.txt", 1536, 1536)