import gzip
import os
import urllib.request
from tkinter import messagebox
from tkinter import ttk


class Downloader():
    # ceck if a file exists and if not asks the user to download it. gets used to check if db are all present
    # TODO: also add checksum check to make sure the index is completely downloaded, if not remove file and download again
    def check_file_exists(self, filepath, url, tk):
        if os.path.exists(filepath):
            return
        else:
            response = messagebox.askquestion("Centrifuge index not found",
                                              "Centrifuge index not found. Do you want to download it?")
            if response == "yes":
                with urllib.request.urlopen(url) as response:
                    # get file size from content-length header
                    file_size = int(response.info().get("Content-Length"))
                    # create progress bar widget
                    progress = ttk.Progressbar(tk, orient="horizontal", length=250, mode="determinate")
                    progress.pack()
                    progress["maximum"] = file_size
                    progress["value"] = 0

                    def update_progress(count, block_size, total_size):
                        progress["value"] = count * block_size
                        tk.update_idletasks()

                    # download file and update progress bar
                    urllib.request.urlretrieve(url, filepath, reporthook=update_progress)
                    progress.destroy()
                messagebox.showinfo("Download complete", "Download complete")

    def unzip_gz(self, file):
        with gzip.open(file, 'rb') as f:
            file_content = f.read()
            file_content.gzip.decompress()
