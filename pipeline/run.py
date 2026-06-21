from rich.console import Console
from rich.panel import Panel
from rich.prompt import IntPrompt

from pipeline.cli.prompts import prompt_create_fcgr_img, prompt_show_fcgr_image

console = Console()

MENU = """
[bold][1][/bold] Create FCGR images
[bold][2][/bold] Show FCGR image
[bold][0][/bold] Exit
"""


def main() -> None:
    console.print(Panel("[bold cyan]FCGR-CNN Pipeline[/bold cyan]", expand=False))

    while True:
        console.print(Panel(MENU.strip(), title="Menu", expand=False))
        choice = IntPrompt.ask("Choice", default=0)

        if choice == 0:
            break
        elif choice == 1:
            prompt_create_fcgr_img()
        elif choice == 2:
            prompt_show_fcgr_image()
        else:
            console.print("[red]Invalid choice. Enter 0, 1, or 2.[/red]")


if __name__ == "__main__":
    main()
