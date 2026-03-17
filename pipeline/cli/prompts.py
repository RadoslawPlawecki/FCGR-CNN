from __future__ import annotations

from rich.console import Console
from rich.panel import Panel
from rich.prompt import Confirm, IntPrompt, Prompt
from rich.table import Table

from pipeline.image_creation.fcgr_utils import create_fcgr_img, show_fcgr_image

console = Console()

VALID_SIZES = {"500", "1000", "1500", "2000"}


def parse_stop_after(value: str) -> int | None:
    if value.strip() == "":
        return None
    n = int(value.strip())  # raises ValueError if not an integer
    if n <= 0:
        raise ValueError(f"stop_after must be a positive integer, got {n}")
    return n


def is_valid_kmer(k: int) -> bool:
    return 2 <= k <= 6


def prompt_create_fcgr_img() -> None:
    console.print(Panel("[bold]Create FCGR Images[/bold]", expand=False))

    while True:
        k_mer = IntPrompt.ask("k-mer length", default=6)
        if is_valid_kmer(k_mer):
            break
        console.print("[red]k-mer must be between 2 and 6.[/red]")

    while True:
        raw = Prompt.ask("Stop after N rows [blank = all]", default="")
        try:
            stop_after = parse_stop_after(raw)
            break
        except ValueError as e:
            console.print(f"[red]{e}[/red]")

    save_original = Confirm.ask("Save ref and mut images?", default=False)

    table = Table(show_header=False, box=None)
    table.add_row("k_mer", str(k_mer))
    table.add_row("stop_after", str(stop_after) if stop_after is not None else "all")
    table.add_row("save_original", str(save_original))
    console.print(Panel(table, title="Parameters", expand=False))

    if not Confirm.ask("Run?", default=True):
        return

    try:
        create_fcgr_img(k_mer=k_mer, stop_after=stop_after, save_original=save_original)
    except Exception:
        console.print_exception()


def prompt_show_fcgr_image() -> None:
    console.print(Panel("[bold]Show FCGR Image[/bold]", expand=False))

    while True:
        k_mer = IntPrompt.ask("k-mer length")
        if is_valid_kmer(k_mer):
            break
        console.print("[red]k-mer must be between 2 and 6.[/red]")

    chromosome = Prompt.ask("Chromosome (accession)")
    location = Prompt.ask("Location")
    label = Prompt.ask("Label")
    size = Prompt.ask("Window size", choices=["500", "1000", "1500", "2000"])

    table = Table(show_header=False, box=None)
    table.add_row("k_mer", str(k_mer))
    table.add_row("chromosome", chromosome)
    table.add_row("location", location)
    table.add_row("label", label)
    table.add_row("size", size)
    console.print(Panel(table, title="Parameters", expand=False))

    if not Confirm.ask("Run?", default=True):
        return

    show_fcgr_image(k_mer=k_mer, chromosome=chromosome, location=location, label=label, size=size)
