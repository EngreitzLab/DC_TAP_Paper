import marimo

__generated_with = "0.16.5"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # DC-TAP-seq Element Category Analysis Notebook

    This notebook generates profile plots for significant WTC11 DE-G pairs stratified by


    This notebook explores the significant WTC11 DE-G pairs where the assigned element 
    category is "CTCF Element".


    Targeted Questions: 

    Questions:
    1. Do they (CTCF sites) show some low level of H3K27ac?
    2. Are the guides that exhibit the strongest effects at these
       elements very close to the CTCF binding site (motif)?
    """
    )
    return


@app.cell
def _():
    import os
    return


if __name__ == "__main__":
    app.run()
