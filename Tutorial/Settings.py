class Settings:
    def __init__(self):
        self.directory='2a'
        self.use_vcd=True

        self.shift_factor=0.98
        self.cutoff_absolute=False
        self.cutoff=0.015
        self.sigma_1=0.192
        self.sigma_2=0.1
        self.x_min=1000
        self.x_max=1500
        self.w=12
    def get(self):
        return self.shift_factor, self.cutoff_absolute, self.cutoff, self.sigma_1, self.sigma_2, self.use_vcd, self.x_min, self.x_max
    def get_directory(self):
        return self.directory
    def get_w(self):
        return self.w


