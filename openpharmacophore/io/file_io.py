import abc


class AbstractIO(abc.ABC):

    def __init__(self, filename: str):
        self.filename = filename

    def __enter__(self) -> "AbstractIO":
        return self

    def __exit__(self, *args):
        pass

    @abc.abstractmethod
    def load_data(self):
        raise NotImplementedError
