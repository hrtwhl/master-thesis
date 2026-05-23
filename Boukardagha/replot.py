import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

# Sicherstellen, dass das aktuelle Verzeichnis im Pfad ist, um config/data zu finden
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    from config import ASSET_NAMES, ensure_dirs, OUTPUT_DIR
    from data import load_all
    from reporting import export_all, _apply_style
    from backtest import BacktestResult
    print(">>> Module erfolgreich importiert.")
except ImportError as e:
    print(f">>> FEHLER beim Import: {e}")
    sys.exit(1)

def replot_only():
    print(">>> Starte replot_only Funktion...")
    
    # 1. Verzeichnisse prüfen
    ensure_dirs()
    print(f">>> Zielverzeichnis: {OUTPUT_DIR}")

    # 2. Basisdaten laden
    print(">>> Lade Basisdaten via load_all()... (dies kann kurz dauern)")
    try:
        bundle = load_all()
        print(f">>> Basisdaten geladen. Returns Shape: {bundle['returns'].shape}")
    except Exception as e:
        print(f">>> FEHLER beim Laden der Basisdaten: {e}")
        return

    # 3. CSV Pfad prüfen und laden
    csv_path = Path(OUTPUT_DIR) / "tables" / "daily_backtest_output.csv"
    print(f">>> Suche CSV unter: {csv_path}")
    
    if not csv_path.exists():
        print(f">>> FEHLER: Datei {csv_path} nicht gefunden!")
        return

    try:
        df = pd.read_csv(csv_path, index_col="date", parse_dates=True)
        print(f">>> CSV erfolgreich geladen ({len(df)} Zeilen).")
    except Exception as e:
        print(f">>> FEHLER beim Lesen der CSV: {e}")
        return

    # 4. Result-Objekt rekonstruieren
    print(">>> Rekonstruiere BacktestResult-Objekt...")
    try:
        # Extrahiere Gewichtsspalten (w_SPX, w_BOND, etc.)
        weights_cols = [c for c in df.columns if c.startswith("w_")]
        weights_df = df[weights_cols].rename(columns=lambda x: x[2:])
        
        # Validierung der Assets
        missing_assets = [a for a in ASSET_NAMES if a not in weights_df.columns]
        if missing_assets:
            print(f">>> WARNUNG: Assets fehlen in CSV: {missing_assets}")

        reconstructed_result = BacktestResult(
            pnl=df["pnl"],
            weights=weights_df,
            cum_pnl=df["cum_pnl"],
            K_history=pd.Series(df["K"], index=df.index, name="K"),
            tpl_label=pd.Series(df["regime"], index=df.index, name="regime"),
            tpl_max_prob=pd.Series(df["max_p"], index=df.index, name="max_p"),
            tpl_count=pd.Series(df["G"], index=df.index, name="G"),
            turnover=df["turnover"]
        )
        print(">>> BacktestResult-Objekt erfolgreich erstellt.")
    except KeyError as e:
        print(f">>> FEHLER: Erwartete Spalte in CSV fehlt: {e}")
        return

    # 5. Styling und Export
    print(">>> Wende Styles an (_apply_style)...")
    _apply_style()
    
    print(">>> Starte export_all (Generierung der PNGs)...")
    try:
        artefacts = export_all(reconstructed_result, bundle["returns"])
        print(f">>> ERFOLG: {len(artefacts['figures'])} Grafiken neu erstellt.")
    except Exception as e:
        print(f">>> FEHLER während export_all: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    print(">>> Skript manuell gestartet.")
    replot_only()
    print(">>> Skript beendet.")