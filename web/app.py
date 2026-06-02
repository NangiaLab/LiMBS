from flask import Flask, render_template, request
from lipid_registry import LIPIDS
from cg_templates import generate_cg_block

app = Flask(__name__)


@app.route("/")
def home():
    return render_template("form.html", lipids=LIPIDS)


@app.route("/generate", methods=["POST"])
def generate():
    mem_type = request.form.get("type", "planar").strip().lower()
    protein = request.form.get("protein", "").strip()

    box_x = request.form.get("box_x", "42").strip()
    box_y = request.form.get("box_y", "42").strip()
    box_z = request.form.get("box_z", "15").strip()

    sol = request.form.get("sol", "W").strip()
    salt_conc = request.form.get("salt_conc", "0.15").strip()
    salt_ion = request.form.get("salt_ion", "NaCl").strip()
    rand = request.form.get("rand", "0.1").strip()
    solr = request.form.get("solr", "0.5").strip()

    selected_lipids = []
    upper_parts = []
    lower_parts = []

    for lipid in LIPIDS:
        upper = request.form.get(f"{lipid}_upper", "0").strip()
        lower = request.form.get(f"{lipid}_lower", "0").strip()

        try:
            upper_val = float(upper)
            lower_val = float(lower)
        except ValueError:
            upper_val = 0.0
            lower_val = 0.0

        if upper_val > 0 or lower_val > 0:
            selected_lipids.append(lipid)

        if upper_val > 0:
            upper_parts.append(f"{lipid}:{upper}")

        if lower_val > 0:
            lower_parts.append(f"{lipid}:{lower}")

    if not selected_lipids:
        return """
        <html>
        <body style="font-family: Arial, sans-serif; margin: 40px;">
            <h2>No lipids selected</h2>
            <p>Please enter at least one non-zero lipid value.</p>
            <p><a href="/">Go back</a></p>
        </body>
        </html>
        """

    cg_blocks = " || ".join(
        generate_cg_block(
            lipid,
            LIPIDS[lipid]["class"],
            LIPIDS[lipid]["tail_a"],
            LIPIDS[lipid]["tail_b"],
        )
        for lipid in selected_lipids
    )

    if mem_type == "vesicle":
        headers = "; ".join(
            f"{LIPIDS[lipid]['header']}, APL:{LIPIDS[lipid].get('apl', 0.65)}"
            for lipid in selected_lipids
        )

        vesicle_upper_block = ", ".join(f"{lipid}:1" for lipid in selected_lipids)
        vesicle_lower_block = ", ".join(f"{lipid}:1" for lipid in selected_lipids)

        notation = (
            f"{headers} || "
            f"leaflets -u{{{vesicle_upper_block}}} -l{{{vesicle_lower_block}}} || "
            f"type:vesicle || "
            f"box:[{box_x},{box_y},{box_z}] || "
            f"{cg_blocks}"
        )

        next_step = """
        <b>Vesicle system detected</b>

        <p>
        If you are using Martini 3 for a vesicle system, use the
        LiMBS-to-TS2CG parser script. The parser reads the generated
        LiMBS line notation and converts it into the corresponding
        structure and topology files.
        </p>
        """

    else:
        headers = "; ".join(LIPIDS[lipid]["header"] for lipid in selected_lipids)

        upper_block = ", ".join(upper_parts)
        lower_block = ", ".join(lower_parts)

        protein_part = f"; Protein:{protein}" if protein else ""

        notation = (
            f"{headers}{protein_part} || "
            f"leaflets -u{{{upper_block}}} -l{{{lower_block}}} || "
            f"type:planar || "
            f"box:[{box_x},{box_y},{box_z}] || "
            f"sol:{sol} || "
            f"salt:{salt_conc} {salt_ion} || "
            f"rand:{rand} || "
            f"solr:{solr} || "
            f"{cg_blocks}"
        )

        next_step = """
        <b>Planar membrane detected</b>

        <p>
        If you are using Martini 3 for a planar membrane, use the
        LiMBS-to-INSANE parser script. The parser reads the generated
        LiMBS line notation and converts it into the corresponding
        structure and topology files.
        </p>
        """

    return f"""
    <html>
    <head>
    <style>
    body {{
        font-family: Arial, sans-serif;
        margin: 40px;
        line-height: 1.5;
    }}

    h1, h2, h3 {{
        color: #222222;
    }}

    .notation-box {{
        background-color: #f5f5f5;
        border: 1px solid #cccccc;
        border-radius: 8px;
        padding: 15px;
        margin-top: 10px;
        margin-bottom: 30px;
        overflow-x: auto;
    }}

    .workflow-box {{
        background-color: #fafafa;
        border-left: 4px solid #4a90e2;
        padding: 15px;
        margin-top: 10px;
        margin-bottom: 30px;
    }}

    .note-box {{
        background-color: #fffde7;
        border-left: 4px solid #f4c542;
        padding: 15px;
        margin-top: 10px;
        margin-bottom: 25px;
    }}

    pre {{
        white-space: pre-wrap;
        word-wrap: break-word;
        font-size: 14px;
    }}

    a {{
        color: #1a5fb4;
        text-decoration: none;
    }}

    a:hover {{
        text-decoration: underline;
    }}
    </style>
    </head>

    <body>

    <h1>Generated LiMBS Notation</h1>

    <div class="notation-box">
        <pre>{notation}</pre>
    </div>

    <h2>Recommended Next Step</h2>

    <div class="workflow-box">
        {next_step}
    </div>

    <div class="note-box">
        <b>Note:</b><br><br>
        This prototype generates LiMBS line notation only.<br><br>
        External LiMBS parser scripts can be used separately to convert the
        generated notation into simulation-ready structure and topology files.<br><br>
        Parser execution is not currently integrated into the web interface.
    </div>

    <p>
        <a href="/">Generate another LiMBS notation</a>
    </p>

    </body>
    </html>
    """


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
