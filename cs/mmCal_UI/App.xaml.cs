using MmCalculatorWpf;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Windows;

namespace mmCal_UI {
    public partial class App : Application {
        protected override void OnStartup(StartupEventArgs e) {
            try {
                base.OnStartup(e);

                // 埋め込みDLLを展開してロード
                string dllPath = NativeDllLoader.LoadEmbeddedDll();

                // MainWindow を明示生成（StartupUriは使わない方針）
                var win = new MainWindow(dllPath); // DLL パスを渡す
                MainWindow = win;
                win.Show();
            } catch (Exception ex) {
                MessageBox.Show(ex.ToString(), "Startup Error");
                Shutdown();
            }
        }
    }

    internal static class NativeDllLoader {
        public static string? DllPath { get; private set; }
        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern IntPtr LoadLibrary(string lpFileName);

        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern bool SetDllDirectory(string lpPathName);

        public static string LoadEmbeddedDll() {
            var asm = Assembly.GetExecutingAssembly();

            // -----------------------------
            // プロセスのアーキテクチャ判定
            // -----------------------------
            bool is64 = Environment.Is64BitProcess;
            string archFolder = is64 ? "x64" : "x86";
            string dllName = is64 ? "mmCal_x64.dll" : "mmCal_x86.dll";
            string resourceName = "mmCal_UI." + dllName;

            // -----------------------------
            // 展開先フォルダ
            // -----------------------------
            string tempDir = Path.Combine(Path.GetTempPath(), "mmCalculator", archFolder);
            Directory.CreateDirectory(tempDir);
            string dllPath = Path.Combine(tempDir, "mmCal_UI.dll"); // 固定名で統一

            // -----------------------------
            // バージョンチェック
            // -----------------------------
            bool shouldExtract = true;

            if (File.Exists(dllPath)) {
                var existingVersion = FileVersionInfo.GetVersionInfo(dllPath).FileVersion;

                string tempCheck = Path.Combine(tempDir, "check_mmCal_UI.dll");
                using (Stream? s = asm.GetManifestResourceStream(resourceName)) {
                    if (s == null) throw new Exception("Embedded DLL not found: " + resourceName);
                    using (var fs = new FileStream(tempCheck, FileMode.Create, FileAccess.Write))
                        s.CopyTo(fs);
                }

                var embeddedVersion = FileVersionInfo.GetVersionInfo(tempCheck).FileVersion;
                File.Delete(tempCheck);

                if (existingVersion == embeddedVersion) shouldExtract = false;
                else {
                    try { File.Delete(dllPath); } catch { } // 他プロセスでロック中なら無視
                }
            }

            // -----------------------------
            // 埋め込みDLLを展開
            // -----------------------------
            if (shouldExtract) {
                using Stream? s = asm.GetManifestResourceStream(resourceName);
                if (s == null) throw new Exception("Embedded DLL not found: " + resourceName);

                using var fs = new FileStream(dllPath, FileMode.Create, FileAccess.Write);
                s.CopyTo(fs);
            }

            // -----------------------------
            // DLL フォルダを検索パスに追加
            // -----------------------------
            if (!SetDllDirectory(tempDir)) {
                int err = Marshal.GetLastWin32Error();
                throw new Exception($"SetDllDirectory failed with error {err}");
            }

            // -----------------------------
            // DLLロード（安全のため残す）
            // -----------------------------
            IntPtr h = LoadLibrary(dllPath);
            if (h == IntPtr.Zero) {
                int err = Marshal.GetLastWin32Error();
                MessageBox.Show($"LoadLibrary failed: {err}\nDLL Path: {dllPath}", "DLL Load Error", MessageBoxButton.OK, MessageBoxImage.Error);
                throw new Exception("LoadLibrary failed: " + err);
            } else {
                //MessageBox.Show($"DLL loaded successfully!\nHandle: 0x{h:X}\nPath: {dllPath}", "DLL Load OK");
            }
            DllPath = dllPath;
            return dllPath; // DllImport に渡す用
        }
    }
}
