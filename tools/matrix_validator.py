import numpy as np
import re
import time
import datetime

def load_array(path):
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()

    text = text.replace("{", "[")
    text = text.replace("}", "]")

    # I → j
    text = text.replace("I", "j")

    # Float表記整理
    text = re.sub(r"\+ ", "+", text)
    text = re.sub(r"\- ", "-", text)

    # numpy直接パース
    arr = np.fromstring(
        text.replace("[", "").replace("]", ""),
        dtype=np.complex128,
        sep=","
    )

    # 形状復元
    size = int(np.sqrt(len(arr)))
    arr = arr.reshape((size, size))

    return arr


# 誤差評価
def compare_arrays(file1, file2, eps=1e-15):
    now = datetime.datetime.now() 
    today = now.strftime("'%y/%m/%d/%H:%M:%S") 
    print(today+" compare started...")
    start = time.perf_counter()
    a = load_array(file1)
    print("Loaded: " + file1)
    b = load_array(file2)
    end = time.perf_counter()
    elapsed = (end - start)
    print("Loaded: " + file2)
    print(f"Load time: {elapsed:.6f}" + "[s]")
    
    

    if a.shape != b.shape:
        raise ValueError(f"Shape mismatch: {a.shape} vs {b.shape}")

    diff = a - b
    abs_err = np.abs(diff)

    # 基本誤差
    max_err = np.max(abs_err)
    mean_err = np.mean(abs_err)

    # 相対誤差
    denom = np.abs(b) + eps
    rel_err = abs_err / denom

    max_rel_err = np.max(rel_err)
    mean_rel_err = np.mean(rel_err)

    # エネルギー保存誤差
    energy_a = np.sum(np.abs(a)**2)
    energy_b = np.sum(np.abs(b)**2)

    energy_err = np.abs(energy_a - energy_b)

    # 結果表示
    print("===================================")
    print("Shape:", a.shape)
    print("-----------------------------------")

    print("Max Absolute Error :", max_err)
    print("Mean Absolute Error:", mean_err)

    print("Max Relative Error :", max_rel_err)
    print("Mean Relative Error:", mean_rel_err)

    print("Energy A:", energy_a)
    print("Energy B:", energy_b)
    print("Energy Conservation Error:", energy_err)
    
    end = time.perf_counter()
    elapsed = (end - start)
    print(f"\nTotal elapsed time: {elapsed:.6f} [s]")
    print("===================================")

    return {
        "max_err": max_err,
        "mean_err": mean_err,
        "max_rel_err": max_rel_err,
        "mean_rel_err": mean_rel_err,
        "energy_err": energy_err
    }


# -------------------------------------------------------
if __name__ == "__main__":    
    print("===================================")
    print("Matrix Validator  -mmKreutzef")
    print("===================================")
    compare_arrays(
        "FFT_tool\\fft_random_1536.txt",
        "FFT_tool\\out_random.txt"
    )