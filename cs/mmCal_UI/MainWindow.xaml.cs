using System;
using System.Collections.Generic;
using System.Text;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media;
using System.Runtime.InteropServices;

namespace MmCalculatorWpf
{
    public partial class MainWindow : Window
    {
        private const string PlaceholderText = "";
        private int evalCount = 0;

        private List<string> inputHistory = new();
        private List<string> outputHistory = new();

        public MainWindow()
        {
            InitializeComponent();
            InputHistoryListBox.ItemsSource = inputHistory;
            OutputHistoryListBox.ItemsSource = outputHistory;
            ShowPlaceholder();
            Loaded += (s, e) =>
            {
                InputTextBox.Focus();            // WPF フォーカス
                Keyboard.Focus(InputTextBox);    // キーボードフォーカス
            };
        }

        private void ShowPlaceholder()
        {
            InputTextBox.Text = PlaceholderText;
            //InputTextBox.Foreground = Brushes.Gray;
            InputTextBox.Foreground = Brushes.White;
        }

        private void HidePlaceholder()
        {
            if (InputTextBox.Text == PlaceholderText)
            {
                InputTextBox.Text = "";
                InputTextBox.Foreground = Brushes.White;
            }
        }

        private void InputTextBox_TextChanged(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {
            if (string.IsNullOrEmpty(InputTextBox.Text))
                PlaceholderTextBlock.Visibility = Visibility.Visible;
            else
                PlaceholderTextBlock.Visibility = Visibility.Collapsed;
        }

        private void InputTextBox_GotFocus(object sender, RoutedEventArgs e) => HidePlaceholder();
        private void InputTextBox_LostFocus(object sender, RoutedEventArgs e)
        {
            if (string.IsNullOrWhiteSpace(InputTextBox.Text))
                ShowPlaceholder();
        }

        private void InputTextBox_KeyDown(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Enter && Keyboard.Modifiers == ModifierKeys.Shift)
            {
                e.Handled = true; // Shift+Enter で改行しない
                EvaluateInput();
            }
        }

        private void EvaluateInput()
        {
            string expr = InputTextBox.Text.Trim();
            if (string.IsNullOrEmpty(expr) || expr == PlaceholderText)
                return;

            int outNeed, errNeed, errPos;
            var outSb = new StringBuilder(1024);
            var errSb = new StringBuilder(256);

            int status = NativeMethods.mmcal_eval_ex(
                NativeMethods.ctx,
                expr,
                outSb, outSb.Capacity, out outNeed,
                errSb, errSb.Capacity, out errNeed,
                out errPos
            );

            if (status == 0)
                AddHistory(expr, outSb.ToString());
            else
                ShowError(expr, errSb.ToString(), errPos);

            InputTextBox.Clear();
            ShowPlaceholder();
        }

        private void AddHistory(string input, string output)
        {
            evalCount++;
            inputHistory.Add($"In[{evalCount}]: {input}");
            outputHistory.Add($"Out[{evalCount}]: {output}");
            InputHistoryListBox.Items.Refresh();
            OutputHistoryListBox.Items.Refresh();
            RecentOutputTextBox.Text = output;
        }

        private void ShowError(string expr, string error, int pos)
        {
            string pointer = pos >= 0 && pos <= expr.Length ? new string(' ', pos) + "^" : "";
            RecentOutputTextBox.Text = $"Error: {error}\n{expr}\n{pointer}";
        }

        private void EvaluateButton_Click(object sender, RoutedEventArgs e)
        {
            EvaluateInput();
        }

        //private async void HistoryListBox_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        //{
        //    if (sender is System.Windows.Controls.ListBox listBox && listBox.SelectedItems.Count > 0)
        //    {
        //        // 選択アイテムを改行区切りで結合してコピー
        //        var sb = new StringBuilder();
        //        foreach (var item in listBox.SelectedItems)
        //            sb.AppendLine(item.ToString());

        //        Clipboard.SetText(sb.ToString());

        //        // 通知表示
        //        CopyNotificationTextBlock.Text = "Copied to clipboard!";
        //        CopyNotificationTextBlock.Visibility = Visibility.Visible;

        //        // 2秒後に非表示
        //        await Task.Delay(750);
        //        CopyNotificationTextBlock.Visibility = Visibility.Collapsed;
        //    }
        //}






        private void ClearHistory()
        {
            inputHistory.Clear();
            outputHistory.Clear();
            InputHistoryListBox.Items.Refresh();
            OutputHistoryListBox.Items.Refresh();
            evalCount = 0;
            RecentOutputTextBox.Clear();
        }
    }

    internal static class NativeMethods
    {
        [DllImport("mmCal_x64.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern IntPtr mmcal_create();

        [DllImport("mmCal_x64.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern void mmcal_destroy(IntPtr ctx);

        [DllImport("mmCal_x64.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
        internal static extern int mmcal_eval_ex(
            IntPtr ctx,
            string expr,
            StringBuilder outStr, int outCap, out int outNeed,
            StringBuilder errStr, int errCap, out int errNeed,
            out int errPos);

        [DllImport("mmCal_x64.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
        internal static extern int mmcal_last_error(IntPtr ctx, StringBuilder outStr, int outCap);

        internal static IntPtr ctx = mmcal_create();
    }
}
