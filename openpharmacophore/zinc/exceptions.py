

class ZincDownloadError(IOError):
    """ Exception raised when a file format is not supported or
        is incorrect.
    """

    def __init__(self, status_code, url):

        self.message = f"Error downloading file from {url}.\nStatus code: {status_code}"
        super().__init__(self.message)
