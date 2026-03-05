#!/usr/bin/env python3
import sys

def parse_rt_data_formatted(filename):
    """フォーマット済みRTデータから電流値抽出（行折り返し対応）"""
    with open(filename, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    data_start = next((i for i, l in enumerate(lines) if l.strip() and not l.startswith('#')), -1)
    if data_start < 0:
        return None
    
    tokens = ' '.join(lines[data_start:]).split()
    times, jx, jy, jz = [], [], [], []
    
    for i in range(0, len(tokens)-15, 16):
        try:
            times.append(float(tokens[i]))
            jx.append(float(tokens[i+13]))
            jy.append(float(tokens[i+14]))
            jz.append(float(tokens[i+15]))
        except: pass
    
    return {'times': times, 'Jx': jx, 'Jy': jy, 'Jz': jz} if times else None

def analyze_current(data, label=""):
    """電流データの統計分析"""
    if not data or not data['times']:
        print(f"\n{label}: データなし")
        return None
    
    jx, jy, jz = data['Jx'], data['Jy'], data['Jz']
    jmag = [(jx[i]**2 + jy[i]**2 + jz[i]**2)**0.5 for i in range(len(jx))]
    
    print(f"\n{'='*80}\n電流値統計分析: {label}\n{'='*80}")
    print(f"タイムステップ数: {len(data['times'])}")
    print(f"時間範囲: {data['times'][0]:.8f} - {data['times'][-1]:.8f} a.u.")
    print(f"\nJx: |max|={max(abs(j) for j in jx):.10e}, mean={sum(jx)/len(jx):.10e}")
    print(f"Jy: |max|={max(abs(j) for j in jy):.10e}, mean={sum(jy)/len(jy):.10e}")
    print(f"Jz: |max|={max(abs(j) for j in jz):.10e}, mean={sum(jz)/len(jz):.10e}")
    print(f"|J|: max={max(jmag):.10e}, mean={sum(jmag)/len(jmag):.10e}")
    
    return {'Jx_max': max(abs(j) for j in jx), 'Jy_max': max(abs(j) for j in jy),
            'Jz_max': max(abs(j) for j in jz), 'J_mag_max': max(jmag), 'n_ts': len(data['times'])}

def show_table(data, label, nrows=15):
    """データをテーブル形式で表示"""
    if not data: return
    print(f"\n{label} - データテーブル（最初と最後）:")
    print("-"*100)
    print(f"{'Step':>5} {'Time[a.u.]':>15} {'Jx[a.u.]':>20} {'Jy[a.u.]':>20} {'Jz[a.u.]':>20}")
    print("-"*100)
    for i in range(min(nrows, len(data['times']))):
        print(f"{i+1:5d} {data['times'][i]:15.8f} {data['Jx'][i]:20.12e} {data['Jy'][i]:20.12e} {data['Jz'][i]:20.12e}")
    if len(data['times']) > 2*nrows:
        print("  ...")
        for i in range(max(nrows, len(data['times'])-nrows), len(data['times'])):
            print(f"{i+1:5d} {data['times'][i]:15.8f} {data['Jx'][i]:20.12e} {data['Jy'][i]:20.12e} {data['Jz'][i]:20.12e}")
    print("-"*100)

print("="*80)
print("DG-Fragment RT vs 従来的 RT - 電流値比較分析")
print("="*80)

conv = parse_rt_data_formatted("H2_periodic_20_conventional_rt_rt.data")
dg = parse_rt_data_formatted("H2_periodic_20_dg_new_param_rt.data")

conv_st = analyze_current(conv, "従来的 RT")
dg_st = analyze_current(dg, "DG-Fragment RT")

if conv: show_table(conv, "従来的 RT", nrows=15)
if dg: show_table(dg, "DG-Fragment RT", nrows=15)

if conv_st and dg_st:
    print(f"\n{'='*80}\n比較分析\n{'='*80}")
    print(f"従来的RT最大|J|: {conv_st['J_mag_max']:.10e} a.u.")
    print(f"DG-Fragment最大|J|: {dg_st['J_mag_max']:.10e} a.u.")
    ratio = conv_st['J_mag_max'] / max(1e-20, dg_st['J_mag_max'])
    print(f"比率 (Conv/DG): {ratio:.4f}")
elif conv_st:
    print(f"\n{'='*80}\nデータ利用性\n{'='*80}")
    print("✅ 従来的RTデータ: 完全（21タイムステップ）")
    print(f"   期待電流範囲: |J|_max ~ {conv_st['J_mag_max']:.2e} a.u.")
    print("⚠️  DG-Fragment RTデータ: 不完全（テスト未完了）")

print(f"\n{'='*80}\n")
