#!/usr/bin/env python3
"""Render a standalone HTML review from the synthetic lot report."""
import argparse
import json
from html import escape
from pathlib import Path


def esc(value):
    return escape(str(value), quote=True)


def render(report):
    qc_labels = {
        'expected_positive_detected': 'Témoin positif attendu',
        'negative_control_clear': 'Témoin négatif',
        'minimum_depth_met': 'Profondeur des échantillons',
    }
    cards = ''.join(
        f'<div class="check {"good" if ok else "bad"}"><span>{esc(qc_labels[key])}</span>'
        f'<strong>{"CONFORME AU SCÉNARIO" if ok else "ÉCHEC"}</strong></div>'
        for key, ok in report['qc'].items()
    )
    sample_rows = ''.join(
        f'<tr><td>{esc(sid)}</td><td>{esc(info["role"])}</td><td>{esc(info["total_reads"]):}</td>'
        f'<td>{esc(info.get("expected_taxon", "—"))}</td></tr>'
        for sid, info in report['samples'].items()
    )
    candidates = ''.join(
        '<tr>'
        f'<td><strong>{esc(w["sample_id"])}</strong><small>{esc(w["taxon"])} · {esc(w["reference_accession"])}</small></td>'
        f'<td>{esc(w["reads"])}</td><td>{esc(w["unique_regions"])}</td>'
        f'<td>{esc(round(w["breadth"]*100,2))} %</td>'
        f'<td>{esc(w["identity_pct"])} %</td><td>{esc(w["host_similarity_pct"])} %</td>'
        f'<td>{esc(", ".join(w["review_flags"]) or "Aucun indicateur complémentaire")}</td>'
        f'<td><span class="pill">{esc(w["state"])}</span></td></tr>'
        for w in report['worklist']
    ) or '<tr><td colspan="8">Aucun candidat au seuil de démonstration.</td></tr>'
    hits = ''.join(
        '<tr>'
        f'<td>{esc(h["sample_id"])}</td><td>{esc(h["role"])}</td>'
        f'<td>{esc(h["taxon"])}<small>{esc(h["reference_accession"])}</small></td>'
        f'<td>{esc(h["reads"])}</td><td>{esc(h["rpm"])}</td><td>{esc(h["unique_regions"])}</td>'
        f'<td>{esc(round(h["breadth"]*100,2))} %</td><td>{esc(h["mean_identity"])} %</td>'
        f'<td>{"OUI" if h["gate_pass"] else "NON"}</td></tr>'
        for h in report['hits']
    )
    p = report['provenance']
    status = report['status']
    heading = 'Témoins acceptés : revue humaine possible' if status == 'REVIEW_READY' else 'Lot bloqué : échec d’un témoin ou de la profondeur'
    return f'''<!doctype html>
<html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(report['run_id'])} · NGS Lot Review</title>
<style>
:root {{font-family:system-ui,-apple-system,Segoe UI,sans-serif;color:#173348;background:#f2f6f8}}
* {{box-sizing:border-box}} body {{margin:0}} main {{max-width:1120px;margin:auto;padding:32px 20px 65px}}
small {{display:block;color:#56707e}} h1 {{font-size:clamp(2rem,4vw,3rem);line-height:1.1;margin:8px 0}} h2 {{font-size:1.25rem;margin:0 0 15px}}
.eyebrow {{text-transform:uppercase;letter-spacing:.12em;font-size:.75rem;font-weight:700;color:#457384}}
.sub {{color:#49677a}} .hero {{margin:24px 0;background:#123e53;color:white;padding:26px;border-radius:16px}}
.hero.block {{background:#754325}} .hero strong {{font-size:1.5rem}} .hero p {{margin-bottom:0;color:#e4eef1}}
.card {{background:white;border:1px solid #dce7ea;border-radius:15px;padding:23px;margin:17px 0;box-shadow:0 3px 12px #1733480a}}
.checks {{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}} .check {{padding:13px;border-radius:9px;display:flex;flex-direction:column;gap:5px}}
.good {{background:#eaf8f1}} .bad {{background:#fff0e6}} .check strong {{font-size:.72rem;letter-spacing:.03em}}
.scroll {{overflow-x:auto}} table {{border-collapse:collapse;width:100%;font-size:.91rem}} th,td {{text-align:left;border-bottom:1px solid #e4ecef;padding:11px;vertical-align:top;white-space:nowrap}}
th {{color:#4e6776;background:#f1f6f7}} .pill {{display:inline-block;border-radius:20px;background:#fff0d9;color:#68410b;padding:3px 8px;font-weight:700;font-size:.76rem}}
.note {{padding:14px 18px;border-left:4px solid #cc8840;background:#fff8e9}} code {{overflow-wrap:anywhere;font-size:.79rem}} .footer {{color:#526b79;font-size:.83rem}}
@media(max-width:720px) {{.checks {{grid-template-columns:1fr}} main {{padding:20px 12px}} .card {{padding:15px}}}}
</style></head><body><main>
<div class="eyebrow">Portefeuille technique · données entièrement synthétiques</div>
<h1>NGS Lot Review</h1>
<p class="sub">Lot {esc(report['run_id'])} · {esc(report['platform'])} · {esc(report['assay_context'])}</p>
<section class="hero {'block' if status == 'QC_BLOCKED' else ''}"><strong>{esc(heading)}</strong>
<p>Statut : {esc(status)} · {len(report['worklist'])} candidat(s) de revue. Aucun résultat ne vaut identification confirmée.</p></section>
<section class="card"><h2>1. Portes de contrôle</h2><div class="checks">{cards}</div></section>
<section class="card"><h2>2. Échantillons du lot</h2><div class="scroll"><table><thead><tr><th>ID</th><th>Rôle</th><th>Reads totaux</th><th>Cible attendue</th></tr></thead><tbody>{sample_rows}</tbody></table></div></section>
<section class="card"><h2>3. Liste de revue humaine</h2>
<p class="sub">Ces signaux dépassent les seuils illustratifs. Les indicateurs supplémentaires demandent une expertise, sans éliminer automatiquement le signal.</p>
<div class="scroll"><table><thead><tr><th>Échantillon et référence</th><th>Reads</th><th>Régions</th><th>Couverture</th><th>Identité</th><th>Similarité hôte</th><th>Indicateurs</th><th>État</th></tr></thead><tbody>{candidates}</tbody></table></div></section>
<section class="card"><h2>4. Toutes les preuves fournies</h2><div class="scroll"><table><thead><tr><th>Échantillon</th><th>Rôle</th><th>Taxon et référence</th><th>Reads</th><th>Reads/million</th><th>Régions</th><th>Couverture</th><th>Identité</th><th>Seuil</th></tr></thead><tbody>{hits}</tbody></table></div>
<details><summary>Paramètres de démonstration</summary><pre>{esc(json.dumps(report['demo_rules'], ensure_ascii=False, indent=2))}</pre></details></section>
<section class="card"><h2>5. Traçabilité et limites</h2><p class="note">{esc(report['limitations'])}</p>
<p>Références : <code>{esc(p['reference_snapshot'])}</code> · Logiciel : <code>{esc(p['software'])}</code></p>
<details><summary>Empreintes SHA-256 des fichiers utilisés</summary><p>Hits : <code>{esc(p['hits_sha256'])}</code><br>Métadonnées : <code>{esc(p['metadata_sha256'])}</code><br>Règles : <code>{esc(p['config_sha256'])}</code></p></details></section>
<p class="footer">Prototype indépendant ; aucune donnée client, procédure propriétaire ou méthode iDTECT® utilisée.</p>
</main></body></html>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('results/report.json'))
    parser.add_argument('--output', type=Path, default=Path('results/report.html'))
    args = parser.parse_args()
    report = json.loads(args.input.read_text(encoding='utf-8'))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render(report), encoding='utf-8')
    print(args.output)


if __name__ == '__main__':
    main()
