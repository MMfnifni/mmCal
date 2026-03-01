import math

def write_fft_sine_matrix(
    filename,
    rows,
    cols,
    f1=4, f2=8, f3=2, f4=3,
    A1=1.0, A2=0.7, A3=0.5
):
    with open(filename, "w", encoding="utf-8") as f:
        f.write("{")

        for i in range(rows):
            f.write("{")
            for j in range(cols):

                x = i / rows
                y = j / cols

                val = (
                    A1 * math.sin(2 * math.pi * f1 * x) +
                    A2 * math.sin(2 * math.pi * f2 * y) +
                    A3 * math.sin(2 * math.pi * (f3 * x + f4 * y))
                )

                f.write(f"{val:.10f}")
                if j != cols - 1:
                    f.write(",")

            f.write("}")
            if i != rows - 1:
                f.write(",")

        f.write("})")

    print(f"Written: {filename}")


write_fft_sine_matrix("fft_sin_16.txt", 16, 16)