import tkinter as tk
from tkinter import ttk
import threading
import queue

class ProgressWindow:
    def __init__(self, sample_name):
        self.sample_name = sample_name
        self.root = None
        self.status_label = None
        self.progress_bar = None
        self.stop_button = None
        self.stop_flag = threading.Event()
        self.message_queue = queue.Queue()
        self.running = False
        
    def create_window(self):
        """Create the window in the main thread"""
        if not self.root:
            self.root = tk.Toplevel()
            self.root.title(f"Analysis Progress - {self.sample_name}")
            self.root.geometry("500x200")
            
            # Center the window on screen
            screen_width = self.root.winfo_screenwidth()
            screen_height = self.root.winfo_screenheight()
            x = (screen_width - 500) // 2
            y = (screen_height - 200) // 2
            self.root.geometry(f"500x200+{x}+{y}")
            
            # Make window stay on top
            self.root.attributes('-topmost', True)
            
            # Create main frame with padding
            main_frame = ttk.Frame(self.root, padding="20")
            main_frame.pack(fill='both', expand=True)
            
            # Sample name header
            self.header_label = ttk.Label(
                main_frame, 
                text=f"Analyzing sample: {self.sample_name}",
                font=("Helvetica", 12, "bold")
            )
            self.header_label.pack(pady=(0, 10))
            
            # Current step label
            self.status_label = ttk.Label(
                main_frame, 
                text="Initializing...",
                wraplength=460,
                justify='center'
            )
            self.status_label.pack(pady=(0, 15))
            
            # Progress bar
            self.progress_bar = ttk.Progressbar(
                main_frame,
                mode='indeterminate',
                length=400
            )
            self.progress_bar.pack(pady=(0, 15))
            
            # Stop button
            self.stop_button = ttk.Button(
                main_frame,
                text="Stop Analysis",
                command=self.stop_analysis,
                style='Accent.TButton'
            )
            self.stop_button.pack(pady=(0, 10))
            
            # Configure style for the stop button
            style = ttk.Style()
            style.configure('Accent.TButton', 
                          font=('Helvetica', 10),
                          padding=5)
            
            # Handle window close
            self.root.protocol("WM_DELETE_WINDOW", self.on_close)
            
            self.running = True
            self.progress_bar.start(10)
            self.update_status()
            
    def start(self):
        """Start the progress window in the main thread"""
        if not self.running:
            # Schedule window creation in main thread
            if hasattr(self.root, 'after'):
                self.root.after(0, self.create_window)
            else:
                # If we don't have a root yet, we're probably in the main thread
                self.create_window()
        
    def stop(self):
        """Stop and close the progress window"""
        if self.running:
            self.running = False
            if self.root:
                self.root.after(0, self.root.destroy)
                self.root = None
        
    def update_status(self):
        """Update the status label with new messages from the queue"""
        if not self.running or not self.root:
            return
            
        try:
            while True:
                message = self.message_queue.get_nowait()
                if self.status_label:
                    self.status_label.config(text=message)
        except queue.Empty:
            pass
        
        if self.running and self.root:
            self.root.after(100, self.update_status)
            
    def set_status(self, message):
        """Set a new status message"""
        if self.running:
            self.message_queue.put(message)
            # If we're in the main thread and have a root, update immediately
            if self.root and threading.current_thread() is threading.main_thread():
                self.update_status()
        
    def stop_analysis(self):
        """Set the stop flag and update status"""
        self.stop_flag.set()
        self.set_status("Stopping analysis... Please wait...")
        if self.stop_button:
            self.stop_button.config(state='disabled')
        
    def on_close(self):
        """Handle window close button"""
        self.stop_analysis()
        
    def is_stopped(self):
        """Check if analysis should be stopped"""
        return self.stop_flag.is_set()
