# BOLOS Python 3 実行ガイド

> 公開用フォルダでは、断面積と BOLSIG+ 出力をリポジトリ直下の `Phelps/` に配置します。付属の実行用ノートブックは `run_basic.ipynb` です。原版ソースと比較スクリプト `validate_port.py` は公開対象に含めません。過去の検証結果は参考記録として残しています。文書中の `run_basic.py` は、ご自身で保存して使うコード例です。

このガイドは、提供された Phelps の N₂・O₂ 断面積を使い、電子エネルギー分布関数（EEDF）、反応速度係数、輸送係数を計算するための説明です。

作業場所は、クローンしたリポジトリの直下です。`pyproject.toml`、`bolos/`、`Phelps/` が同じフォルダにある場所でコマンドを実行します。Windows、macOS、Linux で手順は共通で、違うのは仮想環境の作成と有効化だけです。

まず既存のサンプルまたはノートブックを動かし、その後に「5. 基本的な実行プログラム」を読むと、各処理の役割を確認できます。

## 目次

1. [ファイルと実行方法の関係](#1-ファイルと実行方法の関係)
2. [Python 実行環境の確認](#2-python-実行環境の確認)
3. [既存サンプルを実行する](#3-既存サンプルを実行する)
4. [ノートブックを実行する](#4-ノートブックを実行する)
5. [基本的な実行プログラム](#5-基本的な実行プログラム)
6. [プログラムの各部分の意味](#6-プログラムの各部分の意味)
7. [コマンドライン引数の仕組み](#7-コマンドライン引数の仕組み)
8. [結果の保存と電場の変更](#8-結果の保存と電場の変更)
9. [計算結果を確認するポイント](#9-計算結果を確認するポイント)
10. [エラー・警告への対処](#10-エラー警告への対処)
11. [検証の再実行と関連資料](#11-検証の再実行と関連資料)

## 1. ファイルと実行方法の関係

主な構成は次のとおりです。

```text
Py3bolos/                               リポジトリ直下
├─ .venv/                               ローカルに作る Python 仮想環境
├─ bolos/                               計算ライブラリ本体
├─ Phelps/
│  ├─ Phelps_N2-O2_20261026.txt         入力する断面積
│  └─ N2-O2_80-20_2000Td.dat            BOLSIG+ 計算結果
├─ samples/
│  ├─ single.py                         1条件を計算するサンプル
│  └─ bolsig.py                         数値表形式の BOLSIG+ 出力との比較
├─ run_basic.ipynb                      実行用ノートブック
└─ EXECUTION_GUIDE_js.md                このガイド
```

`bolos/solver.py` は計算ライブラリです。通常は、このファイルを直接実行するのではなく、`single.py` などの実行プログラムから読み込みます。

| 実行方法 | 入力・条件の与え方 | 用途 |
|---|---|---|
| `samples/single.py` | コマンドの末尾にファイル名、`--en`、`--temp` を指定 | 1条件をすぐに計算する |
| `run_basic.ipynb` | セル内の変数を変更 | 計算過程・図・表を順に確認する |
| 自分で書く実行プログラム | Python 変数またはコマンドライン引数 | 条件走査や結果保存を組み込む |

`Phelps_N2-O2_20261026.txt` は衝突断面積です。`N2-O2_80-20_2000Td.dat` は BOLSIG+ の計算済み出力で、断面積入力には使用しません。

## 2. Python 実行環境の確認

以下のコマンドは、リポジトリ直下のターミナルで実行します。Python の `>>>` プロンプトに入力するコードとは区別してください。Python は 3.10 以上を使います。macOS と Linux では `python3`、Windows では `py -3` が、用途に合う実行ファイルです。

### 2.1 仮想環境を作って有効にする

macOS と Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Windows（PowerShell）:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
```

有効化できると、プロンプト先頭に `(.venv)` が表示されます。使用する Python を確認します。

```bash
python -c "import sys; print(sys.executable)"
```

出力が、このリポジトリ内の `.venv` を指していれば、対象の仮想環境を使っています。macOS と Linux では `.venv/bin/python`、Windows では `.venv\Scripts\python.exe` です。

Windows で `Activate.ps1` が実行ポリシーにより拒否された場合は、有効化を省略し、以降の `python` を次の実行ファイルに置き換えます。

```powershell
.\.venv\Scripts\python.exe -c "import sys; print(sys.executable)"
```

### 2.2 BOLOS をインストールして読み込めるか確認する

有効化した仮想環境で、リポジトリ直下からインストールします。

```bash
python -m pip install -U pip
python -m pip install -e ".[notebook]"
python -c "from bolos import solver; print(solver.__file__)"
```

`-e` は、このフォルダを参照する形でのインストールです。`[notebook]` は Matplotlib と ipykernel も追加する指定です。引用符は、シェルが `.[notebook]` をファイル名の展開として扱わないために付けます。`python -m pip` と書くことで、実行に使う Python に対してインストールできます。

`solver.__file__` の出力が、このリポジトリ内の `bolos/solver.py` を指していることを確認します。描画だけのサンプルなら `".[samples]"`、BOLSIG+ 比較ノートブックなら `".[validation]"`、計算だけなら追加指定なしの `python -m pip install -e .` を使います。

## 3. 既存サンプルを実行する

仮想環境を有効にしたあと、以降のコマンドはリポジトリ直下で実行します。パス区切りは `/` です。Windows の PowerShell でも同じ表記で実行できます。

### 3.1 そのまま実行できるコマンド

```bash
python samples/single.py Phelps/Phelps_N2-O2_20261026.txt --en=100 --temp=300
```

有効化していない Windows では、`python` の代わりに `.\.venv\Scripts\python.exe` を使います。

各部分の意味は次のとおりです。

| 指定 | 意味 |
|---|---|
| `python` | プログラムを動かす Python |
| `single.py` | 実行するプログラム |
| `Phelps_N2-O2_20261026.txt` のパス | 必須引数 `input`。読み込む断面積ファイル |
| `--en=100` | 換算電場 E/N = 100 Td |
| `--temp=300` | ガス温度 T = 300 K |

`input` という文字そのものを入力する必要はありません。その位置に実際のファイルのパスを書きます。実行後にファイル名を入力する方式ではなく、起動時に指定する方式です。

このコマンドは、現在のターミナルの作業フォルダがリポジトリ直下であることを前提にしています。別のフォルダで実行するときは、パスも合わせて変更してください。空白のあるパスは引用符で囲みます。

### 3.2 条件を変える

200 Td・300 K で計算する例です。

```bash
python samples/single.py Phelps/Phelps_N2-O2_20261026.txt --en=200 --temp=300
```

`single.py` の混合比はコード中で N₂ = 0.8、O₂ = 0.2 に固定されています。混合比を変える場合は該当行を編集するか、ノートブックを使用してください。

引数の一覧と説明は次のコマンドで表示できます。

```bash
python samples/single.py --help
```

`--debug` を付けると、反復計算などの詳細ログが表示されます。

```bash
python samples/single.py Phelps/Phelps_N2-O2_20261026.txt --en=100 --temp=300 --debug
```

### 3.3 エディタの実行ボタンとの関係

実行ボタンから起動したコマンドが `python single.py` だけの場合、断面積ファイルの引数が渡されません。その場合は次のエラーになります。

```text
the following arguments are required: input
```

引数を必要とする `single.py` は、まず上記のコマンドをターミナルに貼り付けて実行してください。「5」のコード例は、条件をプログラム内で設定するため、ファイル名の引数を付けずに実行できます。

### 3.4 `bolsig.py` を使う場合

`samples/bolsig.py` は比較用で、`input` と `bolsigdata` の2つの必須引数を受け取ります。

このサンプルは `np.loadtxt()` で読み込める数値表を想定し、実際のコードでは第2列を電場、第4列以降を各過程の速度係数として参照しています。

提供された `N2-O2_80-20_2000Td.dat` は説明行・条件・EEDF などを含むブロック形式なので、このサンプルに直接渡す形式とは異なります。このブロック形式を扱う原版比較スクリプトはローカル専用で、公開対象から除外しています。

## 4. ノートブックを実行する

1. リポジトリ直下の `run_basic.ipynb` を開きます。
2. カーネルに、BOLOS をインストールした Python 環境を選択します。
3. ノートブックの作業ディレクトリをリポジトリ直下に設定します。`Path.cwd()` が `bolos/` と `Phelps/` のあるフォルダを指すことを確認してください。
4. 最初のセルから上から順番に実行します。

ターミナルの仮想環境と、ノートブックのカーネルは別に選択されます。入力パスは次の形です。

```python
script_dir = Path.cwd()
input_file = script_dir / "Phelps" / "Phelps_N2-O2_20261026.txt"
```

付属ノートブックでは、主に次の変数を変更します。

```python
en_td_values = np.linspace(200.0, 1000.0, 9)
gas_temperature_k = 300.0
fractions = {"N2": 0.6, "O2": 0.4}
```

混合比の合計は 1 にします。条件を変更した後は、そのセル以降を再実行してください。以下の文書にある単一条件のコード例や過去の検証記録とは設定が異なる場合があります。

## 5. 基本的な実行プログラム

BOLOS を呼び出すプログラムは、次の順番で書きます。

```text
ライブラリの読み込み
    ↓
断面積ファイル・計算条件の設定
    ↓
エネルギー格子とソルバの作成
    ↓
断面積の読み込み
    ↓
混合比・温度・電場を設定して init()
    ↓
初期 EEDF を作り converge()
    ↓
必要に応じて格子を作り直して再計算
    ↓
平均エネルギー・輸送係数・反応速度係数の出力
```

以下は、コードの中に条件を書く完全な例です。使用する場合は、**リポジトリ直下に `run_basic.py` という名前で新規保存**してください。このガイドのコード例であり、`run_basic.py` 自体はガイド作成時には追加していません。

この例では、サンプル `single.py` から独立して、ノートブックと同じ格子設定を用いています。

```python
from pathlib import Path

import numpy as np
import scipy.constants as co
import matplotlib.pyplot as plt

from bolos import parser, solver, grid


def main():
    # 1. 入力ファイルと計算条件
    # このプログラムをリポジトリ直下に保存した場合のパスです。
    script_dir = Path(__file__).resolve().parent
    input_file = script_dir / "Phelps" / "Phelps_N2-O2_20261026.txt"

    en_td = 100.0
    gas_temperature_k = 300.0
    fractions = {"N2": 0.8, "O2": 0.2}

    if not input_file.is_file():
        raise FileNotFoundError(input_file)
    if any(value < 0 for value in fractions.values()) or not np.isclose(sum(fractions.values()), 1.0):
        raise ValueError("混合比は非負とし、合計を1にしてください。")

    # 2. エネルギー格子とソルバを作る
    coarse_grid = grid.LinearGrid(0.0, 60.0, 500)
    bsolver = solver.BoltzmannSolver(coarse_grid)

    # 3. 断面積を読み込む
    with input_file.open(encoding="utf-8") as fp:
        collisions = parser.parse(fp)
    bsolver.load_collisions(collisions)

    # 4. 混合比・ガス温度・換算電場を設定する
    for species, fraction in fractions.items():
        bsolver.set_density(species, fraction)
    bsolver.kT = gas_temperature_k * co.k / co.eV
    bsolver.EN = en_td * solver.TOWNSEND
    bsolver.init()

    # 5. 初期分布から1回目の収束計算
    f_initial = bsolver.maxwell(2.0)
    f_coarse = bsolver.converge(f_initial, maxn=1000, rtol=1e-5)

    # 6. 平均エネルギーに合わせて格子を作り直す
    coarse_mean = bsolver.mean_energy(f_coarse)
    refined_grid = grid.QuadraticGrid(0.0, 15.0 * coarse_mean, 400)
    bsolver.grid = refined_grid
    bsolver.init()

    # 元の格子上の分布を、新しい格子上に補間する
    f_guess = refined_grid.interpolate(f_coarse, coarse_grid)
    eedf, iterations, error = bsolver.converge(
        f_guess, maxn=1000, rtol=1e-6, full=True
    )

    # 7. 結果を確認・出力する
    if not np.isfinite(eedf).all() or np.min(eedf) < -1e-12:
        raise ValueError("EEDF に非有限値または負値があります。")

    print(f"E/N = {en_td:g} Td, T = {gas_temperature_k:g} K")
    print(f"反復回数 = {iterations}, 最終誤差 = {error:.3e}")
    print(f"平均エネルギー = {bsolver.mean_energy(eedf):.8g} eV")
    print(f"移動度 × N = {bsolver.mobility(eedf):.8e} 1/(m V s)")
    print(f"拡散係数 × N = {bsolver.diffusion(eedf):.8e} 1/(m s)")

    print("\n反応速度係数:")
    for target, process in bsolver.iter_all():
        rate = bsolver.rate(eedf, process)
        print(f"{str(process):45s} {rate:.8e} m^3/s")

    # 8. EEDF を描画する
    positive = eedf > 0
    plt.semilogy(bsolver.grid.c[positive], eedf[positive])
    plt.xlabel("Electron energy (eV)")
    plt.ylabel("EEDF (eV$^{-3/2}$)")
    plt.title(f"N2/O2, E/N = {en_td:g} Td")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
```

保存後、リポジトリ直下から実行します。

```bash
python run_basic.py
```

この例では、入力ファイルをコード中で指定しているため、コマンドに断面積ファイル名を付ける必要はありません。

`__file__` は実行プログラム自身の場所を表します。ファイルの場所を基準に入力パスを組み立てると、ターミナルの作業ディレクトリが変わった場合にも対応しやすくなります。

## 6. プログラムの各部分の意味

### 6.1 `import` と `main()`

```python
from bolos import parser, solver, grid
```

- `parser`：BOLSIG+ 形式の断面積ファイルを読み込みます。
- `grid`：エネルギー方向の計算格子を作ります。
- `solver`：EEDF の求解と速度係数・輸送係数の計算を行います。

処理を `main()` という関数にまとめると、実行の入口を明確にできます。

```python
if __name__ == "__main__":
    main()
```

この部分は、そのファイルを `python run_basic.py` として直接実行したときに `main()` を呼びます。別のプログラムから import しただけなら、ここから計算は始まりません。

### 6.2 エネルギー格子

```python
coarse_grid = grid.LinearGrid(0.0, 60.0, 500)
```

0～60 eV を500セルに分割します。セルの境界は `coarse_grid.b`、中心は `coarse_grid.c` です。EEDF の配列はセル中心に対応します。

500セルは、この例の初期設定です。必要な上限エネルギーやセル数は条件に依存します。

### 6.3 断面積の読み込み

```python
with input_file.open(encoding="utf-8") as fp:
    collisions = parser.parse(fp)
bsolver.load_collisions(collisions)
```

`parser.parse()` は衝突過程の辞書のリストを返し、`load_collisions()` がそれを計算に使える過程・ターゲットとして登録します。

今回の断面積ファイルからは43過程が読み込まれます。ファイル中のターゲット名 `N2`、`O2` と、混合比設定で使用する名前は一致させてください。

### 6.4 混合比の設定

```python
bsolver.set_density("N2", 0.8)
bsolver.set_density("O2", 0.2)
```

この API の `density` は、今回の使用方法では各成分のモル分率です。`0.8` は80%を意味します。粒子数密度 [m⁻³] の値を直接入力する場所ではありません。

以下も同じ設定です。

```python
bsolver.target["N2"].density = 0.8
bsolver.target["O2"].density = 0.2
```

### 6.5 温度と電場の単位変換

```python
bsolver.kT = gas_temperature_k * co.k / co.eV
bsolver.EN = en_td * solver.TOWNSEND
```

| 変数 | 意味・単位 |
|---|---|
| `gas_temperature_k` | 入力するガス温度 [K] |
| `bsolver.kT` | ガス温度をエネルギーに換算した値 [eV] |
| `en_td` | 入力する換算電場 [Td] |
| `solver.TOWNSEND` | `1e-21`。Td から V m² への換算係数 |
| `bsolver.EN` | 内部に設定する換算電場 [V m²] |

例えば、100 Td を指定するときは `bsolver.EN = 100.0 * solver.TOWNSEND` です。

### 6.6 `init()` を呼ぶタイミング

```python
bsolver.init()
```

設定した混合比・温度・電場・格子に基づき、求解に必要な係数などを準備します。

これらの条件を変更したときは、求解前に再び `init()` を呼びます。格子を `bsolver.grid = refined_grid` で差し替えた後にも必要です。

### 6.7 初期分布と反復計算

```python
f_initial = bsolver.maxwell(2.0)
```

`2.0` は初期 Maxwell 分布の温度パラメータ kTₑ = 2 eV です。ガス温度の300 Kとは別の値であり、最終的な電子温度や平均エネルギーを2 eVに固定する指定でもありません。

```python
eedf, iterations, error = bsolver.converge(
    f_guess, maxn=1000, rtol=1e-6, full=True
)
```

| 引数・戻り値 | 意味 |
|---|---|
| `f_guess` | 反復計算の初期分布 |
| `maxn` | 反復回数の上限 |
| `rtol` | 隣接する反復分布の差に対する収束しきい値 |
| `full=True` | EEDF とともに反復回数・最終誤差も返す指定 |
| `eedf` | 求めた EEDF の配列 |
| `iterations` | 実際の反復回数 |
| `error` | 最終反復で評価した分布の差 |

この実装の誤差指標は `_norm(abs(f_old - f_new))` です。`rtol` を小さくしても、格子誤差や断面積の不確かさまで同じ割合で小さくなるわけではありません。

収束しなかった場合は `ConvergenceError` が発生します。

### 6.8 格子を作り直す理由

最初は仮のエネルギー範囲で解き、求めた平均エネルギーを使って次の格子を作ります。

```python
refined_grid = grid.QuadraticGrid(0.0, 15.0 * coarse_mean, 400)
```

この例では上限を最初の平均エネルギーの15倍とし、400セルの二次格子で再計算します。二次格子は、この条件では低エネルギー側を細かく扱えます。

格子を変えた場合、以前の EEDF を新しい格子に補間します。

```python
f_guess = refined_grid.interpolate(f_coarse, coarse_grid)
```

単にセル数だけを変えて古い EEDF の配列をそのまま渡すと、格子と分布の対応が崩れます。

### 6.9 出力値の意味

| 呼び出し | 結果 | 単位 |
|---|---|---|
| `bsolver.mean_energy(eedf)` | 平均電子エネルギー | eV |
| `bsolver.mobility(eedf)` | 換算移動度 μN | 1/(m V s) |
| `bsolver.diffusion(eedf)` | 換算拡散係数 DN | 1/(m s) |
| `bsolver.rate(eedf, process)` | 指定過程の速度係数 k | m³/s |

移動度と拡散係数は N を掛けた値です。粒子数密度 N が別途わかっている場合、通常の μ と D はそれぞれ N で割って求めます。

速度係数 `k` は、反応回数そのものではありません。`rate(..., weighted=False)` が既定で、ターゲットのモル分率は掛けません。`weighted=True` は設定したモル分率を掛ける指定であり、実際の粒子数密度まで掛ける処理ではありません。

特定の過程は名前でも指定できます。

```python
k_ionization = bsolver.rate(eedf, "N2 -> N2^+")
```

同じターゲット・生成物名に複数の過程がある場合、文字列検索は既定で最初の一致を使います。全過程を確実に扱うときは、`iter_all()` から得た過程オブジェクトを渡してください。

`iter_all()` は内部で構成した ELASTIC も含む過程を返します。ノートブックの反応速度表は読み込んだ入力過程のリストを使用するため、EFFECTIVE の表示などが異なります。

## 7. コマンドライン引数の仕組み

`single.py` は `argparse` で入力を受け取ります。概念は次のコードです。

```python
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="断面積ファイル")
argparser.add_argument("--en", type=float, default=100.0)
argparser.add_argument("--temp", type=float, default=300.0)
args = argparser.parse_args()

input_file = args.input
en_td = args.en
gas_temperature_k = args.temp
```

このコードが次の入力を受け取るとします。

```text
python single.py cs.txt --en=200 --temp=350
```

対応は次のようになります。

```python
args.input == "cs.txt"
args.en == 200.0
args.temp == 350.0
```

`"input"` のように `--` のない引数は位置で渡す引数です。上の設定では必須です。`--en` と `--temp` は省略可能で、省略時は `default` の値を使います。

自分のプログラムでこの形式を使う場合は、「5」の `input_file`、`en_td`、`gas_temperature_k` に値を代入している部分を、引数から受け取る処理に置き換えます。実際の組み込み方は `samples/single.py` の `parse_args()` と `main()` を参照してください。

## 8. 結果の保存と電場の変更

### 8.1 EEDF をテキストで保存する

「5」の `main()` 内で、EEDF を計算した後に追加する例です。

```python
output_dir = script_dir / "results"
output_dir.mkdir(exist_ok=True)

np.savetxt(
    output_dir / "eedf.txt",
    np.column_stack((bsolver.grid.c, eedf)),
    header="energy_eV eedf_eV_minus_3_over_2",
)
```

1列目がエネルギー [eV]、2列目が EEDF [eV⁻³ᐟ²] です。同じ名前のファイルがある場合は上書きします。条件ごとに残したい場合は、例えば `f"eedf_{en_td:g}Td.txt"` のように名前に電場を含めます。

### 8.2 図を保存する

描画命令の後、`plt.show()` より前に追加します。

```python
plt.savefig(output_dir / f"eedf_{en_td:g}Td.png", dpi=200)
```

この例の `output_dir` は8.1のコードで作成します。

### 8.3 複数の電場で計算する

同じ格子・混合比・温度で電場を変更する基本形は次のとおりです。ここでは、すでに断面積・混合比・温度を設定済みの `bsolver` を使います。

```python
for en_td in [50.0, 100.0, 200.0]:
    bsolver.EN = en_td * solver.TOWNSEND
    bsolver.init()
    initial = bsolver.maxwell(2.0)
    eedf = bsolver.converge(initial, maxn=1000, rtol=1e-6)
    print(en_td, bsolver.mean_energy(eedf))
```

これは電場をループで変える書き方の例です。各条件で格子の上限や解像度が十分かは別に確認します。広い範囲の電場を走査する場合は、各条件で「5」の2段階の格子計算を行う方法も使えます。

## 9. 計算結果を確認するポイント

### 9.1 出力の目安

`run_phelps.ipynb` の保存済み実行結果では、N₂:O₂=80:20、300 K、100 Td に対して次の値でした。

| 量 | 保存済み結果 |
|---|---:|
| 平均電子エネルギー | 2.6272569 eV |
| μN | 1.16865409 × 10²⁴ 1/(m V s) |
| DN | 2.56912809 × 10²⁴ 1/(m s) |

これらはノートブックの格子・収束条件に対する結果です。`single.py` は原版の設定を保ち、2段目が「平均エネルギーの150倍まで・200セル」です。一方、ノートブックとこのガイドの完全例は「15倍まで・400セル」で、許容誤差も異なります。そのため、両者の結果には差が生じます。

### 9.2 収束と格子依存性

反復が収束したことと、格子誤差が十分小さいことは別です。精度を調べる際は、次の変更に対して対象の出力値がどれだけ変わるか比較します。

- 最終格子のセル数を増やす。
- 上限エネルギーを変える。
- 収束しきい値を小さくする。
- 初期分布を変えて、同じ条件の結果を比較する。

EEDF の高エネルギー側は小さい値でも、しきい値の高い反応の速度係数に影響します。平均エネルギーだけでなく、必要な反応の速度係数も比較してください。

### 9.3 EEDF の検査

```python
print("有限値のみ:", np.isfinite(eedf).all())
print("最小値:", np.min(eedf))
print("BOLOS内部の規格化値:", bsolver._norm(eedf))
```

`_norm()` は内部メソッドですが、今回の移植版での検証に使えます。値はおおむね1です。`np.sum(eedf)` だけでは、エネルギーの重みや格子幅を反映した規格化の確認になりません。

## 10. エラー・警告への対処

### `ModuleNotFoundError: No module named 'bolos'`

まず、実際の Python を確認します。

```bash
python -c "import sys; print(sys.executable)"
```

出力がリポジトリ内の `.venv` でなければ、「2.1」の仮想環境を有効にします。対象環境でも import できない場合は、リポジトリ直下でその Python を使ってインストールします。

```bash
python -m pip install -e ".[notebook]"
```

### `the following arguments are required: input`

`single.py` の後ろに断面積ファイルのパスが必要です。「3.1」のコマンドを実行してください。

### `the following arguments are required: input, bolsigdata`

`bolsig.py` を引数なしで実行しています。1条件の計算なら `single.py`、提供された BOLSIG+ データとの比較ならノートブックを使います。

### `FileNotFoundError`

指定したファイルの場所と、現在の作業フォルダを確認します。

macOS と Linux:

```bash
pwd
test -f Phelps/Phelps_N2-O2_20261026.txt && echo yes
```

Windows（PowerShell）:

```powershell
Get-Location
Test-Path .\Phelps\Phelps_N2-O2_20261026.txt
```

ファイルが無いと表示された場合は、作業フォルダをリポジトリ直下にするか、指定パスを修正します。

### `KeyError: 'N2'` または `KeyError: 'O2'`

読み込んだ断面積に、その名前のターゲットがあるか確認します。大文字・小文字や空白を含め、データ中の名前に合わせます。

```python
print(list(bsolver.target))
```

### `ConvergenceError`

最大反復回数までに収束判定を満たしていません。電場・温度の単位、混合比、格子と初期分布の対応、条件変更後の `init()` を確認します。そのうえで、格子範囲・解像度・初期分布・反復回数を調整します。

`maxn` を増やすだけで必ず解決するとは限りません。

### `After substracting INELASTIC from EFFECTIVE ... negative cross-section`

今回の N₂ データでは、原版の EFFECTIVE→ELASTIC 変換時に負値が生じ、ゼロにする処理が通ります。移植版でもその挙動を維持しているため、警告が表示されます。これは Python の構文エラーではありません。

データ処理上の警告であり、計算結果の解釈に関わります。詳細は移植・検証資料も参照してください。

### `divide by zero`、`invalid value` などの `RuntimeWarning`

2026-10-27の修正では、ゼロ流束境界と移動度積分の未使用端点の除算を省略しました。
この部分は結果に使わない演算なので無害です。分布の勾配がゼロ付近の衝突積分は、
0/0をゼロに置き換えるのではなく、安定した局所座標の積分で計算します。

`EEDF tail slopes undefined` は別扱いです。原版のゼロ寄与への置換を維持し、
発生エネルギー範囲と該当する離散電子数割合を表示します。
電子数割合が小さくても高しきい値反応に影響し得るため、格子幅・上限の依存性を
確認するまでは無条件に無視できません。初期格子・最終格子でそれぞれ表示されます。

`Energy boundary ... (safe)` や `is outside ... (safe)` は無害と判断できる情報です。
表示するには `import logging; logging.getLogger().setLevel(logging.INFO)` を実行してください。
断面積クリップの影響範囲は隣接点間の補間も考慮し、今回のN₂では1000～10000 eVです。
新実装でなお RuntimeWarning が出た場合は、原因未確認として扱ってください。

警告だけで計算失敗と決めるのではなく、収束、最終 EEDF の有限性・非負性、規格化、出力係数を確認してください。最終出力に `NaN` や `inf` が残る場合、その結果は使用できません。

### 図が出ない、またはプログラムが終了しない

`plt.show()` は実行環境によって図のウィンドウを表示し、閉じるまで待ちます。図を閉じるとプログラムが終了する場合があります。

画像保存だけを行う場合は、`plt.savefig()` で保存してから `plt.close()` とする方法があります。ノートブックでは通常、図はセルの出力に表示されます。

## 11. 検証の再実行と関連資料

`VALIDATION.md` と `validation_results.json` は過去の検証記録です。原版を必要とするローカル比較スクリプト `validate_port.py` は `.gitignore` で除外されており、原版ソースも公開リポジトリには含めていません。そのため、公開ファイルだけで原版との比較を再実行する構成ではありません。

過去の検証では、Phelps の43過程、N₂:O₂=80:20、300 K で約1・10・100・1000・2000 Td の5条件を比較しました。原版の Python 2 固有の実行操作を互換化し、Python 3 上で計算する方式です。実際の Python 2 環境との直接比較ではありません。通常の計算に原版ソースは不要です。

関連ファイル：

- [付属の基本計算ノートブック](run_basic.ipynb)
- [既存の1条件計算サンプル](samples/single.py)
- [移植内容と制約](PORTING_NOTES.md)
- [過去の検証結果の説明](VALIDATION.md)
- [過去の詳細な検証結果](validation_results.json)
