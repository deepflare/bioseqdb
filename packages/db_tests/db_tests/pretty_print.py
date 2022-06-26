from rich.console import Console

def pretty_print_to_str(*args, **kwargs):
    console = Console()
    with console.capture() as capture:
        console.print(*args, **kwargs)
    str_output = capture.get()
    return str_output[:-1]