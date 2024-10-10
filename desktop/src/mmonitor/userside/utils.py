import tkinter as tk
import customtkinter as ctk

class ToolTip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tip_window = None
        self.alpha = 0
        self.fade_id = None

    def show_tip(self):
        if self.tip_window or not self.text:
            return
        x = self.widget.winfo_rootx() + self.widget.winfo_width() + 5
        y = self.widget.winfo_rooty() + self.widget.winfo_height() // 2

        self.tip_window = tk.Toplevel(self.widget)
        self.tip_window.wm_overrideredirect(True)
        self.tip_window.wm_geometry(f"+{x}+{y}")
        self.tip_window.attributes('-alpha', 0.0)

        style = ctk.get_appearance_mode()
        bg_color = "#2B2B2B" if style == "Dark" else "#F0F0F0"
        fg_color = "#FFFFFF" if style == "Dark" else "#000000"

        label = tk.Label(self.tip_window, text=self.text, justify=tk.LEFT,
                         background=bg_color, foreground=fg_color,
                         relief=tk.SOLID, borderwidth=1,
                         font=("Helvetica", "10", "normal"))
        label.pack(ipadx=5, ipady=5)

        self.fade_in()

    def fade_in(self):
        self.alpha += 0.1
        self.tip_window.attributes('-alpha', min(self.alpha, 1.0))
        if self.alpha < 1.0:
            self.fade_id = self.widget.after(50, self.fade_in)

    def hide_tip(self):
        if self.tip_window:
            if self.fade_id:
                self.widget.after_cancel(self.fade_id)
            self.tip_window.destroy()
            self.tip_window = None
        self.alpha = 0

def create_tooltip(widget, text):
    tooltip = ToolTip(widget, text)
    widget.bind("<Enter>", lambda e: tooltip.show_tip())
    widget.bind("<Leave>", lambda e: tooltip.hide_tip())