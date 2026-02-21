using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Interop;
using System.Windows.Media;
using System.Windows.Threading;

namespace MmCalculatorWpf {
    public partial class MainWindow : Window {
        private int historyIndex = -1;
        private const string PlaceholderText = "";

        private int evalCount = 0;
        private readonly List<HistoryItem> history = new();

        // Evaluate結果を表示している間は、リアルタイム構文チェックを抑制する
        private bool suppressLiveSyntax = false;

        // DLL context
        private IntPtr ctx = IntPtr.Zero;

        internal static string GetDllVersion(string dllPath) {
            if (!File.Exists(dllPath))
                return "unknown";

            var info = FileVersionInfo.GetVersionInfo(dllPath);
            return info.FileVersion ?? "unknown";
        }

        public MainWindow(string dllPath) {
            InitializeComponent();

            NativeMethods.SetDllPath(dllPath);

            HistoryGrid.ItemsSource = history;

            Loaded += async (_, _) => {
                try {
                    ctx = NativeMethods.mmcal_create();
                } catch (Exception ex) {
                    MessageBox.Show("Failed to load mmCal DLL:\n" + ex);
                    Close();
                    return;
                }

                InputTextBox.Focus();
                Keyboard.Focus(InputTextBox);
                ShowPlaceholder();
                Title = $"mmCalculator-UI v{GetDllVersion(dllPath)}";
            };

            Closed += (_, _) => {
                if (ctx != IntPtr.Zero)
                    NativeMethods.mmcal_destroy(ctx);
            };
        }

        // -------------------------
        // Placeholder
        // -------------------------
        private void ShowPlaceholder() {
            InputTextBox.Text = PlaceholderText;
            InputTextBox.Foreground = Brushes.White;
            PlaceholderTextBlock.Visibility = Visibility.Visible;
        }

        private void HidePlaceholder() {
            if (InputTextBox.Text == PlaceholderText) {
                InputTextBox.Text = "";
                InputTextBox.Foreground = Brushes.White;
            }
            PlaceholderTextBlock.Visibility = Visibility.Collapsed;
        }

        private void InputTextBox_GotFocus(object sender, RoutedEventArgs e) => HidePlaceholder();
        private void InputTextBox_LostFocus(object sender, RoutedEventArgs e) {
            if (string.IsNullOrWhiteSpace(InputTextBox.Text))
                ShowPlaceholder();
        }

        private void RecallHistory(int direction) {
            if (history.Count == 0) return;

            // direction: -1=up, +1=down
            if (historyIndex == -1 && direction == -1)
                historyIndex = history.Count - 1;
            else
                historyIndex = Math.Clamp(historyIndex + direction, 0, history.Count - 1);

            InputTextBox.Text = history[historyIndex].Input;
            InputTextBox.CaretIndex = InputTextBox.Text.Length;

            HidePlaceholder(); // プレースホルダーを隠す
        }

        // -------------------------
        // Shift+Enter = Evaluate
        // -------------------------
        private void InputTextBox_KeyDown(object sender, KeyEventArgs e) {
            if (e.Key == Key.Enter && Keyboard.Modifiers == ModifierKeys.Shift) {
                e.Handled = true;
                EvaluateInput();
                historyIndex = -1; // 評価後は履歴ポインタをリセット
                return;
            }

            if (e.Key == Key.Up) {
                e.Handled = true;
                RecallHistory(-1);
            } else if (e.Key == Key.Down) {
                e.Handled = true;
                if (historyIndex == -1) return; // 最新行に戻る
                if (historyIndex < history.Count - 1)
                    RecallHistory(+1);
                else {
                    historyIndex = -1;
                    InputTextBox.Text = "";
                    ShowPlaceholder();
                }
            }
        }


        private void EvaluateButton_Click(object sender, RoutedEventArgs e) => EvaluateInput();

        // -------------------------
        // Live syntax check
        // -------------------------
        private void InputTextBox_TextChanged(object sender, TextChangedEventArgs e) {
            if (suppressLiveSyntax || ctx == IntPtr.Zero)
                return;

            if (string.IsNullOrEmpty(InputTextBox.Text)) {
                PlaceholderTextBlock.Visibility = Visibility.Visible;
                RecentOutputTextBox.Text = "";
                return;
            }

            PlaceholderTextBlock.Visibility = Visibility.Collapsed;

            int errPos;
            var errSb = new StringBuilder(256);

            int st = NativeMethods.mmcal_syntax_check(
                ctx,
                InputTextBox.Text,
                errSb, errSb.Capacity,
                out errPos
            );

            if (st == 0 && errPos < 0) {
                RecentOutputTextBox.Text = "";
                return;
            }

            if (errPos >= 0) {
                string pointerLine = new string(' ', errPos) + "^";
                RecentOutputTextBox.Text =
                    $"Pre-check: {errSb}\n{InputTextBox.Text}\n{pointerLine}";
            } else {
                RecentOutputTextBox.Text = $"Pre-check:\n{errSb}";
            }
        }

        // -------------------------
        // Evaluate core
        // -------------------------
        private async void EvaluateInput() {
            if (ctx == IntPtr.Zero) return;

            string expr = InputTextBox.Text.Trim();
            if (string.IsNullOrEmpty(expr) || expr == PlaceholderText)
                return;

            suppressLiveSyntax = true;

            try {
                int outNeed, errNeed, errPos;
                var outSb = new StringBuilder(1024);
                var errSb = new StringBuilder(256);

                int status = NativeMethods.mmcal_eval_ex(
                    ctx,
                    expr,
                    outSb, outSb.Capacity, out outNeed,
                    errSb, errSb.Capacity, out errNeed,
                    out errPos
                );

                if (status == NativeMethods.MMCAL_OK) {
                    string result = outSb.ToString();
                    AddHistory(expr, result);
                    RecentOutputTextBox.Text = result;
                } else if (status == NativeMethods.MMCAL_S_CLEARED) {
                    ClearHistoryUI();
                    RecentOutputTextBox.Text = "";
                    await ShowToast("Cleared");
                    RecentOutputTextBox.Text = "All histories cleared";
                } else if (status == NativeMethods.MMCAL_S_EXIT) {
                    RecentOutputTextBox.Text = "bye...nara";
                    await ShowToast("bye...nara");
                    await Task.Delay(500);
                    Close();
                    return;
                } else {
                    string err = errSb.ToString();
                    ShowError(expr, err, errPos);
                }

                InputTextBox.Clear();
                ShowPlaceholder();
            } catch (Exception ex) {
                RecentOutputTextBox.Text = ex.ToString();
            } finally {
                await Dispatcher.InvokeAsync(() => { suppressLiveSyntax = false; }, DispatcherPriority.Background);
            }
        }

        // -------------------------
        // History
        // -------------------------
        private void AddHistory(string input, string output) {
            evalCount++;
            history.Add(new HistoryItem { Index = evalCount, Input = input, Output = output });
            HistoryGrid.Items.Refresh();
            if (history.Count > 0) HistoryGrid.ScrollIntoView(history[^1]);
        }

        private void ClearHistoryUI() {
            history.Clear();
            HistoryGrid.Items.Refresh();
            evalCount = 0;
        }

        private void ShowError(string expr, string error, int pos) {
            string pointer = (pos >= 0 && pos <= expr.Length) ? new string(' ', pos) + "^" : "";
            RecentOutputTextBox.Text = $"Error: {error}\n{expr}\n{pointer}";
        }

        // -------------------------
        // Toolbar buttons
        // -------------------------
        private void ClearButton_Click(object sender, RoutedEventArgs e) {
            InputTextBox.Text = "Clear[]";
            EvaluateInput();
        }

        private async void CopyInButton_Click(object sender, RoutedEventArgs e) {
            var item = GetCurrentRowItem();
            if (item == null) { await ShowToast("No cell selected"); return; }
            await CopyToClipboardAsync(item.Input, "Copied In");
        }

        private async void CopyOutButton_Click(object sender, RoutedEventArgs e) {
            var item = GetCurrentRowItem();
            if (item == null) { await ShowToast("No cell selected"); return; }
            await CopyToClipboardAsync(item.Output, "Copied Out");
        }

        private async void CopyBothButton_Click(object sender, RoutedEventArgs e) {
            var item = GetCurrentRowItem();
            if (item == null) { await ShowToast("No cell selected"); return; }
            string text = $"In[{item.Index}]: {item.Input} = {item.Output}";
            await CopyToClipboardAsync(text, "Copied Both");
        }

        private HistoryItem? GetCurrentRowItem() {
            if (HistoryGrid.SelectedCells.Count > 0)
                return HistoryGrid.SelectedCells[0].Item as HistoryItem;
            return HistoryGrid.CurrentCell.Item as HistoryItem;
        }

        private async void HistoryGrid_MouseDoubleClick(object sender, MouseButtonEventArgs e) {
            if (HistoryGrid.CurrentCell.Column == null) return;
            if (HistoryGrid.CurrentItem is not HistoryItem item) return;

            int col = HistoryGrid.CurrentCell.Column.DisplayIndex;
            if (col == 1) await CopyToClipboardAsync(item.Input, "Copied In");
            else if (col == 2) await CopyToClipboardAsync(item.Output, "Copied Out");
        }

        private async Task CopyToClipboardAsync(string text, string toast) {
            if (string.IsNullOrEmpty(text)) { await ShowToast("Empty"); return; }
            IntPtr hwnd = new WindowInteropHelper(this).Handle;
            bool ok = Win32Clipboard.TrySetText(text, hwnd);
            if (ok) await ShowToast(toast);
            else await ShowToast("Clipboard busy...");
        }

        private async Task ShowToast(string msg) {
            CopyToastTextBlock.Text = msg;
            CopyToastTextBlock.Visibility = Visibility.Visible;
            await Task.Delay(700);
            CopyToastTextBlock.Visibility = Visibility.Collapsed;
        }

        private async void CopyRecentButton_Click(object sender, RoutedEventArgs e) {
            string text = RecentOutputTextBox.Text;
            if (string.IsNullOrWhiteSpace(text)) { await ShowToast("No recent output"); return; }
            IntPtr hwnd = new WindowInteropHelper(this).Handle;
            bool ok = Win32Clipboard.TrySetText(text, hwnd);
            if (ok) await ShowToast("Copied Recent");
            else await ShowToast("Clipboard busy...");
        }
    }

    // -------------------------
    // History item
    // -------------------------
    public class HistoryItem {
        public int Index { get; set; }
        public string Input { get; set; } = "";
        public string Output { get; set; } = "";
    }

    // -------------------------
    // Native methods
    // -------------------------
    internal static class NativeMethods {
        private static string dllPath = "mmCal_UI.dll";
        private static IntPtr hModule = IntPtr.Zero;

        internal static void SetDllPath(string path) => dllPath = path;

        // Delegate cache
        private static mmcal_create_delegate? createFunc;
        private static mmcal_destroy_delegate? destroyFunc;
        private static mmcal_eval_ex_delegate? evalExFunc;
        private static mmcal_syntax_check_delegate? syntaxCheckFunc;

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        private delegate IntPtr mmcal_create_delegate();
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        private delegate void mmcal_destroy_delegate(IntPtr ctx);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        private delegate int mmcal_eval_ex_delegate(IntPtr ctx, string expr,
            [Out] StringBuilder outBuf, int outCap, out int outNeed,
            [Out] StringBuilder errBuf, int errCap, out int errNeed,
            out int errPos);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        private delegate int mmcal_syntax_check_delegate(IntPtr ctx, string expr,
            [Out] StringBuilder errBuf, int errCap, out int errPos);

        private static void EnsureLoaded() {
            if (hModule != IntPtr.Zero) return;

            hModule = LoadLibrary(dllPath);
            if (hModule == IntPtr.Zero) throw new Exception("DLL not loaded: " + dllPath);

            createFunc = GetDelegate<mmcal_create_delegate>("mmcal_create");
            destroyFunc = GetDelegate<mmcal_destroy_delegate>("mmcal_destroy");
            evalExFunc = GetDelegate<mmcal_eval_ex_delegate>("mmcal_eval_ex");
            syntaxCheckFunc = GetDelegate<mmcal_syntax_check_delegate>("mmcal_syntax_check");
        }

        private static T GetDelegate<T>(string procName) where T : Delegate {
            IntPtr ptr = GetProcAddress(hModule, procName);
            if (ptr == IntPtr.Zero) throw new Exception($"EntryPoint not found: {procName}");
            return Marshal.GetDelegateForFunctionPointer<T>(ptr);
        }

        internal static IntPtr mmcal_create() { EnsureLoaded(); return createFunc!(); }
        internal static void mmcal_destroy(IntPtr ctx) { EnsureLoaded(); destroyFunc!(ctx); }
        internal static int mmcal_eval_ex(IntPtr ctx, string expr,
            StringBuilder outBuf, int outCap, out int outNeed,
            StringBuilder errBuf, int errCap, out int errNeed, out int errPos) {

            EnsureLoaded();
            return evalExFunc!(ctx, expr, outBuf, outCap, out outNeed, errBuf, errCap, out errNeed, out errPos);
        }
        internal static int mmcal_syntax_check(IntPtr ctx, string expr,
            StringBuilder errBuf, int errCap, out int errPos) {

            EnsureLoaded();
            return syntaxCheckFunc!(ctx, expr, errBuf, errCap, out errPos);
        }

        internal const int MMCAL_OK = 0;
        internal const int MMCAL_S_CLEARED = -13;
        internal const int MMCAL_S_EXIT = -14;

        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern IntPtr LoadLibrary(string lpFileName);
        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern IntPtr GetProcAddress(IntPtr hModule, string procName);
    }

    // -------------------------
    // Win32 Clipboard
    // -------------------------
    internal static class Win32Clipboard {
        private const uint CF_UNICODETEXT = 13;
        private const uint GMEM_MOVEABLE = 0x0002;

        [DllImport("user32.dll", SetLastError = true)]
        private static extern bool OpenClipboard(IntPtr hWndNewOwner);
        [DllImport("user32.dll", SetLastError = true)]
        private static extern bool CloseClipboard();
        [DllImport("user32.dll", SetLastError = true)]
        private static extern bool EmptyClipboard();
        [DllImport("user32.dll", SetLastError = true)]
        private static extern IntPtr SetClipboardData(uint uFormat, IntPtr hMem);
        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern IntPtr GlobalAlloc(uint uFlags, UIntPtr dwBytes);
        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern IntPtr GlobalLock(IntPtr hMem);
        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern bool GlobalUnlock(IntPtr hMem);
        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern IntPtr GlobalFree(IntPtr hMem);

        public static bool TrySetText(string text, IntPtr ownerHwnd, int retry = 10, int delayMs = 25) {
            text ??= "";
            for (int i = 0; i < retry; i++) {
                if (!OpenClipboard(ownerHwnd)) { System.Threading.Thread.Sleep(delayMs); continue; }
                try {
                    EmptyClipboard();
                    byte[] bytes = Encoding.Unicode.GetBytes(text + "\0");
                    IntPtr hGlobal = GlobalAlloc(GMEM_MOVEABLE, (UIntPtr)bytes.Length);
                    if (hGlobal == IntPtr.Zero) return false;

                    IntPtr ptr = GlobalLock(hGlobal);
                    if (ptr == IntPtr.Zero) { GlobalFree(hGlobal); return false; }
                    try { Marshal.Copy(bytes, 0, ptr, bytes.Length); } finally { GlobalUnlock(hGlobal); }

                    IntPtr res = SetClipboardData(CF_UNICODETEXT, hGlobal);
                    if (res == IntPtr.Zero) { GlobalFree(hGlobal); return false; }
                    return true;
                } finally { CloseClipboard(); }
            }
            return false;
        }
    }
}
