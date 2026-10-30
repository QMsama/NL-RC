"""TN施肥响应、背景、投饵经验项与出口负荷。

RCC沿用田面TN代理沟渠出口浓度；负荷按系统出口水量乘浓度计算。
"""
import math
from datetime import date
import numpy as np

ALPHA = {"RM": 0.17910631036786312, "RCC": 0.18862475682454152}
KDECAY = {"RM": 0.32084970234924104, "RCC": 0.436}
BACKGROUND = 0.38139121394675596
GAMMA = 1.0217016137162531
FEED_KG_DAY = 3.25
PROTEIN = 0.32
SYSTEM_HA = 0.31
FEED_CONC = GAMMA * FEED_KG_DAY * PROTEIN / 6.25 / SYSTEM_HA

COEF = {
    "basal": dict(a0=-21.469, aN=0.841, apH=3.074, aSOM=0.0),
    "tillering": dict(a0=-101.2, aN=0.54, apH=22.184, aSOM=-11.093),
}


def c0_from_grid(N_kg_ha, ph, som, event):
    c = COEF[event["stage"]]
    raw = c["a0"] + c["aN"] * (N_kg_ha * event["frac"]) + c["apH"] * ph + c["aSOM"] * som
    return np.maximum(0.0, raw)



def regional_concentration(treatment, day, size, events, peaks):
    concentration = np.full(size, BACKGROUND)
    if treatment == 'RCC':
        concentration = concentration + FEED_CONC
    active = [j for j, event in enumerate(events[treatment]) if day >= event['response']]
    if active:
        j = active[-1]
        elapsed = (day-events[treatment][j]['response']).days
        concentration = concentration + peaks[treatment, j]*math.exp(-KDECAY[treatment]*elapsed)
    return concentration


def point_concentration(treatment, day, events, parameters):
    active = sorted([e for e in events if e['treatment'] == treatment and e['response_date'] <= day],
                    key=lambda e:e['response_date'])
    background = parameters['background_feed_parameters']['background_mg_L']
    feed = 0.0
    if treatment == 'RCC':
        feed = (parameters['background_feed_parameters']['gamma']*parameters['feed_kg_day']
                *parameters['protein_fraction']/6.25/parameters['system_area_ha'])
    fertilizer = 0.0
    if active:
        event = active[-1]
        c = COEF[event['stage']]
        raw = (c['a0']+c['aN']*float(event['N_kg_ha'])+c['apH']*float(event['pH'])
               +c['aSOM']*float(event['SOM_percent_assumed']))
        peak = parameters['peak_parameters'][treatment]
        elapsed = (date.fromisoformat(day)-date.fromisoformat(event['response_date'])).days
        fertilizer = peak['alpha']*max(0.0, raw)*math.exp(-peak['k']*elapsed)
    # 点位基准累加顺序：背景→投饵→施肥。
    total = background+feed
    total += fertilizer
    return {'fertilizer_mg_L':fertilizer, 'background_mg_L':background,
            'feed_mg_L':feed, 'total_mg_L':total}


def outlet_load_kg(outlet_m3, concentration_mg_L):
    return outlet_m3*concentration_mg_L*.001


def validate_regional_events(config, start, end):
    """Validate an explicit per-treatment seasonal nitrogen allocation."""
    if not isinstance(config, dict) or set(config) != {'RM', 'RCC'}:
        raise ValueError('Events must define exactly RM and RCC')
    result = {}
    for treatment, rows in config.items():
        if not isinstance(rows, list) or not rows:
            raise ValueError('Each treatment requires fertilizer events')
        parsed = []
        for row in rows:
            try:
                day = date.fromisoformat(row['date'])
                response = date.fromisoformat(row['response'])
                fraction = float(row['frac'])
                stage = row['stage']
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError('Invalid fertilizer event') from exc
            if not start <= day <= response <= end:
                raise ValueError('Event and response must lie within season in order')
            if stage not in COEF or not math.isfinite(fraction) or fraction < 0:
                raise ValueError('Invalid stage or nitrogen fraction')
            if parsed and (day <= parsed[-1]['date'] or response <= parsed[-1]['response']):
                raise ValueError('Event and response dates must be strictly increasing')
            parsed.append(dict(date=day, response=response, stage=stage, frac=fraction))
        if not math.isclose(sum(e['frac'] for e in parsed), 1.0, rel_tol=0, abs_tol=1e-12):
            raise ValueError('Seasonal nitrogen fractions must sum to one')
        result[treatment] = parsed
    return result
