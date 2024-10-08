import customtkinter as ctk

class DatabaseFrame(ctk.CTkFrame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        ctk.CTkLabel(self, text="Database Management").pack(pady=10)
        ctk.CTkButton(self, text="Manage Emu Database", command=self.manage_emu_db).pack(pady=5)
        ctk.CTkButton(self, text="Manage Centrifuge Database", command=self.manage_centrifuge_db).pack(pady=5)
        ctk.CTkButton(self, text="Go to Analysis", command=lambda: controller.show_frame("AnalysisFrame")).pack(pady=5)

    def manage_emu_db(self):
        # Implement Emu database management logic here
        pass

    def manage_centrifuge_db(self):
        # Implement Centrifuge database management logic here
        pass