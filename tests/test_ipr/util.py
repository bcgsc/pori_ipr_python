class QueryMock:
    def __init__(self, return_values) -> None:
        self.return_values = return_values
        self.index = -1

    def __call__(self, *args, **kwargs):
        self.index += 1
        ret_val = self.return_values[self.index] if self.index < len(self.return_values) else []
        return ret_val
