# figures/

Figures used in the top-level `README.md`. Developer-only; not needed to run
DANTE_TIR.

| file | role |
|---|---|
| `make_workflow_pptx.py` | source of the workflow scheme; draws it with native PowerPoint shapes |
| `dante_tir_workflow.pptx` | generated slide |
| `dante_tir_workflow.png` | generated image, embedded in `README.md` |
| `environment.yml` | conda environment for the script (python-pptx) |

## Editing the workflow scheme

The script is the source. Do not edit the `.pptx` by hand: the changes are
lost on the next regeneration.

Layout follows the DANTE_LTR workflow figure: numbered rows, grey caption,
green status on the right, Okabe-Ito colours. Coordinates are pixels of an
1800 x 1050 canvas; constants at the top of the script hold shared x positions
and colours.

The TSD key in row 4 mirrors the superfamily-specific detectors in
`dt_utils.R` (`detect_tir`, `detect_tir_CACTA`, `detect_tir_pif_harbinger`).
Update it when those constraints change.

## Regenerating

```bash
mamba env create -f figures/environment.yml
mamba activate dante_tir_figures
python figures/make_workflow_pptx.py figures/dante_tir_workflow.pptx
soffice --headless --convert-to pdf --outdir /tmp figures/dante_tir_workflow.pptx
pdftoppm -png -singlefile -scale-to-x 1800 -scale-to-y -1 \
    /tmp/dante_tir_workflow.pdf figures/dante_tir_workflow
```

The PNG is rendered with LibreOffice and poppler (`pdftoppm`). Text uses
DejaVu Sans; install it before rendering, otherwise LibreOffice substitutes
a font and label widths change. Commit the script, `.pptx` and `.png`
together.
