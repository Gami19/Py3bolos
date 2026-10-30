# BOLOS Python 3 翻訳版

> 現在の配布物には、この文書で言及する過去の `run_phelps.ipynb` は含まれていません。付属の計算例は `run_basic.ipynb` と `run_basic.py`、検証には `validate_port.py` を使用してください。

隣接する `bolos-master` の Python 2 コードをもとにした独立したコピーです。
元フォルダ・Phelps データは変更していません。元の著者表記と LGPLv2 の
ライセンス指定を維持しています。公開 API、物理モデル、離散化、反復法、
境界条件、断面積の補間方法は元コードを維持しています。
2026-10-27に未使用境界の演算回避と平坦な分布の衝突積分の安定化を追加しました。

## 変更箇所

- パッケージ内の暗黙の相対 import を明示的な相対 import に変更。
- `fp.next()` → `next(fp)`、`iteritems()` → `items()`、`xrange()` → `range()`。
- `(str, unicode)` の判定を Python 3 の `str` に変更。
- サンプル・チュートリアルの `print` を関数形式に変更。
- 文書設定の除算用特殊メソッドを Python 3 形式に変更。
- サンプルの LaTeX 文字列を raw 文字列に変更。
- `setup.py` の Python 対応表記を変更し、ビルド設定 `pyproject.toml` を追加。
- 廃止された `scipy.integrate.simps` の代わりに `bolos/_compat.py` を追加。
  元の既定値 `even='avg'`（偶数個の点では左右の台形補正を平均）を維持。
  現在の `simpson` への単純置換による計算結果の変化を避けています。
  この互換関数は BOLOS が利用する一次元配列に限定しています。

## インストールと実行

作業ルート `Bolos_analysis` で PowerShell から実行します。

```powershell
& .\env_bolos_analysis\Scripts\python.exe -m pip install -e ".\Py3bolos-master[samples]"
& .\env_bolos_analysis\Scripts\python.exe .\Py3bolos-master\validate_port.py
& .\env_bolos_analysis\Scripts\python.exe .\Py3bolos-master\samples\single.py .\Phelps\Phelps_N2-O2_20261026.txt --en=100 --temp=300
```

通常の利用方法は元コードと同じ `from bolos import parser, solver, grid` です。
実行用ノートブックは `run_phelps.ipynb` です。カーネルに
`env_bolos_analysis/Scripts/python.exe` を選択して上から実行してください。
電場・温度・混合比の設定、EEDF の描画、全反応速度の表示、BOLSIG+ 比較、
任意の JSON 保存を含みます。実行には追加で `ipykernel` が必要です。
Python の最低バージョン表記は 3.10 とし、実機での検証環境は結果 JSON に記録します。
すべての Python 3 バージョンでの動作を検証したという意味ではありません。

## 検証内容と限界

`validate_port.py` は提供された Phelps 断面積 43 過程を読み込み、N2:O2=80:20、
300 K で 1、10、100、1000、2000 Td 付近の BOLSIG+ 記録と同じ電場を使用します。
500 セルの線形格子で初期収束後、上限を平均エネルギーの15倍とする400セルの
二次格子で再収束させます。収束許容誤差は最終段で 1e-6 です。

- 旧 Simpson 法を独立した区間ごとの多項式積分と比較。
- EEDF の有限性・非負性・規格化、全反応速度係数・輸送係数を検査。
- 原版のソルバ・反応過程・ターゲットをメモリ上で読み込み、Python 2 固有の
  実行時操作だけを互換化した参照計算と EEDF・全反応速度・輸送係数を比較。
- 提供された BOLSIG+ 出力との平均エネルギー・移動度・拡散係数の相対差を記録。

結果は `validation_results.json` に保存します。Python 2 実環境での直接比較では
ありません。原版と移植版の比較は同じ現在の NumPy・SciPy・物理定数を使うため、
古いライブラリ環境とのビット単位の一致を保証するものではありません。
BOLSIG+ は独立した実装・格子・設定を持つため、差は翻訳誤りの判定に直結しません。

未使用境界の除算は省略し、平坦・ほぼ平坦な分布の衝突積分は局所座標での
8点 Gauss–Legendre 積分を使います。勾配ゼロでは多項式積分に一致します。
非有限な裾の勾配をゼロ寄与にする原版の処理は維持し、該当する電子数の割合と
エネルギー範囲を警告します。電子数割合だけで反応速度への影響は判断できません。
旧版との完全一致は要求せず、輸送係数の相対差1e-7、反応速度係数1e-4
（絶対許容差1e-28 m³/s）、EEDFは相対1e-6・ピーク値の1e-7の絶対差で比較します。
裾の相対誤差だけでは評価せず、独立積分7ケースと新実装の RuntimeWarning 不在も検証します。
例えば非弾性過程がまったくない場合の空配列処理など、原版の制約は残っています。
文書の Sphinx ビルドと BOLSIG+ 比較用プロットサンプルの全実行は検証対象外です。
