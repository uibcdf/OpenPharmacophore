

class InvalidFileFormatError(ValueError):

    def __init__(self, file_format: str):
        self.message = f"File format {file_format} is not supported"
        super().__init__(self.message)
