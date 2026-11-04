import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

# 1. Загружаем JADES из VizieR (берем оттуда только координаты и z > 9)
jades = pd.read_csv('asu.tsv', sep=';', comment='#', low_memory=False)

for col in ['zphot', 'RA_ICRS', 'DE_ICRS']:
    jades[col] = pd.to_numeric(jades[col], errors='coerce')

jades = jades.dropna(subset=['zphot', 'RA_ICRS', 'DE_ICRS'])
jades_z9 = jades[jades['zphot'] > 9.0]
print(f"В JADES найдено кандидатов с z > 9: {len(jades_z9)}")

# 2. Загружаем скачанные файлы CANDELS (Santini 2015) ради физики
candels_mass = Table.read('hlsp_candels_hst_wfc3_goodss_santini_v1_mass-cat.fits').to_pandas()
candels_phys = Table.read('hlsp_candels_hst_wfc3_goodss_santini_v1_physpar-cat.fits').to_pandas()

candels = pd.merge(candels_mass, candels_phys, on='Seq')

# Умный поиск колонок
mass_cols = [c for c in candels.columns if any(sub in c.lower() for sub in ['m_med', 'mcor', 'mass', 'm_star'])]
sfr_cols = [c for c in candels.columns if 'sfr' in c.lower() and '68' not in c.lower() and 'err' not in c.lower()]
ra_col = next((c for c in candels.columns if c.lower() in ['ra', 'radeg']), 'RA')
dec_col = next((c for c in candels.columns if c.lower() in ['dec', 'decdeg']), 'DEC')

mass_col = mass_cols[0]
sfr_col = sfr_cols[0]

candels_clean = candels[(candels[mass_col] > 0) & (candels[sfr_col] > 0)]

# 3. Кросс-матчинг по координатам
coords_jades = SkyCoord(ra=jades_z9['RA_ICRS'].values*u.degree, dec=jades_z9['DE_ICRS'].values*u.degree)
coords_candels = SkyCoord(ra=candels_clean[ra_col].values*u.degree, dec=candels_clean[dec_col].values*u.degree)

idx, d2d, d3d = coords_jades.match_to_catalog_sky(coords_candels)

# ВОЗВРАЩАЕМ ЖЕСТКИЙ РАДИУС (как было раньше)
match_mask = d2d < 0.5 * u.arcsec

# 4. Формируем итоговую надежную выборку
final_jades = jades_z9[match_mask].copy()
matched_candels = candels_clean.iloc[idx[match_mask]].copy()

final_jades['M_star'] = matched_candels[mass_col].values
final_jades['SFR'] = matched_candels[sfr_col].values
final_jades['sSFR'] = final_jades['SFR'] / final_jades['M_star']

print(f"\nКРОСС-МАТЧИНГ УСПЕШЕН.")
print(f"Итоговая надежная выборка: {len(final_jades)} галактик.")