from typing import Protocol, Union, runtime_checkable


@runtime_checkable
class Interface(Protocol):
    def run_task(self):
        pass


class Factory:
    wrappers = {}
    T = Union[Interface]  # Supported input types

    @classmethod
    def get_interface(cls, input_object):
        if isinstance(input_object, Interface):
            return input_object

        # Search for wrapper
        for wrapper_input, wrapper in Factory.wrappers.items():
            if isinstance(input_object, wrapper_input):
                return wrapper(input_object)

        raise TypeError

    @classmethod
    def register(cls, wrapped_class):
        def decorator(wrapper_class):
            # Add wrapper class to dictionary of wrappers
            Factory.wrappers[wrapped_class] = wrapper_class
            Factory.T = Union[tuple([Interface] + list(Factory.wrappers.keys()))]
            return wrapper_class
        return decorator


def execute(a: Factory.T):
    a = Factory.get_interface(a)
    a.run_task()
