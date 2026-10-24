# 妄想空間 — 素人の像をリーマン予想の最前線に写像する

**Riemann Lab: mapping an amateur's geometric intuitions onto the research frontier of the Riemann Hypothesis**

**English version (main): [README.md](README.md)**

2026年10月15日〜24日。数学の訓練を持たない一人の想像力と、AI（Claude）による翻訳・数値実験の往復で作られた、22 枚の図と対応表の記録です。

> **これは証明ではありません。** ここにある図はすべて、既知の定理・既知の予想・既知の未解決点を、数値で再確認したものです。新しい定理は一つも含まれていません。
> この文書の価値は、「比喩で語られた像が、現代数学のどの部品に正確に対応し、どこで対応が切れるか」を一枚の地図にしたことにあります。

---

## 0. 一分で読む要約

- **出発点の像**：コーヒーカップ（トーラス）の表面に素数の布を巻き、回転するロウソク（光源）で照らすと、影がスクリーンに干渉縞を作る。零点はその縞の節。「光と影が釣り合う線」が実部 1/2。
- **確かめたこと**：この像は、ヴェイユの明示公式（素数の帯）、拡大縮小のユニタリ性（1/2 の出所）、モンゴメリー＝オドリズコの対相関（もつれ）、ベリー＝キーティングの地平線（零点の個数）、ベリーの対角近似（ワームホール＝ランプ）、ハーディ＝リトルウッド／ボゴモルニー＝キーティング（素数対のリズム＝プラトー）、コンヌの吸収スペクトル（穴を編む）、セルバーグ跡公式（双曲編み）、指標の直交性（白色光）、そしてヴェイユの曲線上の証明（光円錐＝ホッジ指数定理）に、部品ごとに対応する。
- **未解決の正体**：ヴェイユの正値性。「どんな模様を織っても、重ねた布の影は負にならない」ことを、零点を使わずに素数側から示すこと。曲線の世界では曲面 X×X の光円錐（符号 (+,−,−,…)）がそれを保証した。整数の世界では、その**曲面に当たる空間**（本記録では「妄想空間」と呼ぶ）が未定義。
- **到達点**：部品と設計図は（ヴェイユ 1948 のものとして）すべて像に対応づいた。欠けているのは組み立てではなく、部品を置く場所。

---

## 1. 像 ↔ 数学 の対応表

| 像（コンノ） | 数学の対象 | 状態 | 図 |
|---|---|---|---|
| コーヒーカップ（トーラス） | 非可換トーラス A_θ、モジュラー曲面、E×E | 定義済み | — |
| 素数の帯・輪 | 明示公式の閉軌道項 Λ(n)/√n、長さ log p | **定理**（ヴェイユ） | 1 |
| 回転する光源、常に正面から同じ明るさ | 拡大縮小の流れ、ユニタリ性 ⇒ 指数 1/2 | **定理** | 3 |
| 光と影の均衡線 | 関数等式 ξ(s)=ξ(1−s)、折り目 Re s=1/2 | **定理**（テイト：積公式から） | — |
| 偽の光（軸外の零点） | e^{δT} で成長する項 | **定理**（明示公式） | 2 |
| 螺旋 | 固有振動 x^{iγ}=e^{iγ log x}、メリン変換 | **定理** | 4 |
| もつれ（近くで強く、遠くで弱い） | 対相関 1−(sin πr/πr)² | **予想**（モンゴメリー）、数値一致 | 5, 19 |
| ホログラフィー／ワームホール | 形状因子のランプ＝素数の対角和 | **半古典**（ベリー 1985）、数値一致 | 6, 15 |
| ブラックホール | 地平線つき xp 模型の計数 N(T) | **半古典**（ベリー＝キーティング）、数値一致 | 14 |
| 密度の二乗、カップ二つを重ねる | 自己相関 → ヴェイユの二次形式 | **定理**（RH と同値） | 7 |
| ファンデルワールス力（層の和） | — | **外れ**：素数は和ではない | 8 |
| チェーンステッチ（一本抜くとほどける） | オイラー積（乗法的）／穴は局所的 | **半分**：壊れ方は抜いた素数で決まる | 8, 9 |
| 目を増やして編む | 双曲面、セルバーグ跡公式、重み ℓ/(2 sinh(ℓ/2)) | **定理**（符号が逆） | 11 |
| 裏返して見る／裏編み／穴を編む | 関数等式／奇数次（ドニンガー）／補空間（コンヌ） | 定理／計画／定理 | 13 |
| トーラスの中心を通るワームホール、均衡 | M(x), ψ(x)−x の幅が x^{1/2} | **RH と同値**（向きは言えない：ポリア反例） | 18 |
| 前の素数に素数二つを足すと次の素数 | ゴールドバッハ（隙間版）、表し方の数 = S(h) | **予想**、1000万まで反例なし | 18 |
| 螺旋の隙間に大きな素数が座る | エラトステネスの篩、メルテンス、隙間模様の自己相関 = S(h) | **定理**＋**予想**（HL） | 17 |
| 色つき螺旋、光の重ね合わせ | ディリクレ指標、直交性＝白色光 | **定理** | 19, 20 |
| 垂直に交わる二つのトーラス | E×E 上の繊維 H, V：H²=V²=0, H·V=1（光円錐の母線） | **定理**（曲線） | 22 |
| 影がぴたりと一致する場所 | 対角線 Δ、レフシェッツ不動点公式 | **定理** | 21 |
| 出会いは点だが線になる | 閉軌道（長さ log p）、流れの跡公式（ギルマン、コンヌ） | **定理** | 1, 21 |
| 光円錐の垂線への射影 | Γ'=Γ−pH−V（極の項を引く）、ホッジ指数定理 | **定理**（曲線）／**未解決**（整数） | 21, 22 |
| 1/2 はトーラスの真ん中の糸 | 時間軸 H+V ＝ ゼータの極 s=0,1 | 定理（曲線）＋数値（図 8） | 8, 22 |
| **妄想空間** | ℤ×ℤ の曲面（コンヌ＝コンサーニの算術サイト等が候補） | **未定義** | — |

---

## 2. 22 の実験

各図は `figures/`、スクリプトは `scripts/`、零点データは `data/` にあります。括弧内は検算の精度。

1. **explicit_formula_orbits.png** — ヴェイユの明示公式。零点 300 個だけから、素数の閉軌道（長さ log p^k、重み Λ(n)/√n）が再構成される（差 1e-13）。
2. **fake_zero.png** — 軸上／軸外の偽零点を混ぜる。軸外の零点は e^{δT} で成長し帳簿を壊すが、素数側の背景 e^{T/2} の方が速く成長するため遠い偽零点は素朴には見えない（難しさの正体）。
3. **unitary_half.png** — (U_t f)(x)=e^{σt}f(e^t x) がユニタリになる指数は σ=1/2 のみ（光源が「同じ明るさで正面から」回る条件）。
4. **prime_spiral.png** — 素数ベクトル Λ(n)n^{-1/2}e^{-iγ log n} の和は零点の周波数でだけ共鳴する。
5. **zero_correlation.png** — 零点 1500 個の間隔分布と対相関は GUE に乗る（モンゴメリー＝オドリズコ）。
6. **form_factor.png** — 形状因子にランプ（傾き 0.89）とプラトー（1.01）。
7. **weil_positivity.png** — 帯状模様の「影の行列」を素数だけから構成。固有値はすべて ≥0（ヴェイユ正値性の有限窓版）、軸外の偽零点を混ぜると負になる。
8. **vdw_test.png** — 素数の手を一本抜く／指数を 1/2 からずらす／一律倍率、いずれも正値性が壊れる。素数は和ではなく全体で一つの構造。
9. **chain_direction.png** — 窓 L=2→8 で、手を一本抜いた壊れ方は飽和し、抜いた素数のオイラー因子の反共鳴 2 log p/(√p+1) で決まる（比 0.90〜0.96）。
10. **quiet_band.png** — 明示公式を周波数で分解。零点のない帯域 0<g<14 で無限遠の光源と素数の合計が正確に相殺（2.5e-4）。素数同士が互いの沈みを埋める。
11. **hyperbolic_crochet.png** — 目を毎列一定割合で増やす編み物＝双曲面。セルバーグの輪の重み ℓ/(2 sinh(ℓ/2)) と素数の重み Λ(n)/√n が ℓ=log p で一致。違いは輪の項の符号だけ。
12. **torus_curvature.png** — 楕円断面のトーラス：外側（正曲率）と穴側（負曲率）の総量が +4π/−4π で厳密に釣り合う（ガウス＝ボネ）。曲率の符号は跡公式の符号とは別物。
13. **knit_holes.png** — 穴を編む：模様 f を整数倍 x,2x,3x,… に重ねた布 F=Σf(nx) を臨界線上で分解すると、どんな模様でも同じ周波数で光が消え、それが零点。MF=ζ·Mf（1e-14）。
14. **horizon_count.png** — xp 模型に地平線（切り落とし）を入れた半古典計数が零点 1500 個の個数と一致。残りのゆらぎ S(T) は平均 0、|S|≤1.23。
15. **wormhole_from_primes.png** — ランプを素数だけから対角近似で組み立て：τ<1 で零点の測定と相関 0.989。τ>1 では対角和は伸び続ける。
16. **prime_pairs_plateau.png** — 1000 万までの素数対を数え、ハーディ＝リトルウッド S(h) と比 0.99〜1.00。形状因子を素数から厳密に組み立て（0.08%）、対角／非対角に分解：非対角は τ<1 でゼロ、τ>1 で負になりプラトーを作る。非対角を S(h) だけで置き換えた模型がプラトーを再現。
17. **spiral_gaps.png** — 2,3,5,7 の螺旋を周長 210 のカップに巻く：7 より大きい素数はすべて隙間に座る。隙間の割合はメルテンス e^{−γ}/log p。隙間模様の自己相関 = S(h)。
18. **center_and_dream.png** — M(x)/√x, (ψ(x)−x)/√x が ±0.5〜0.8 に収まる（RH ⇔ 幅が x^{1/2}）。ポリア予想は 906,150,257 で破れる（向きは言えない）。隙間は素数一つか二つの和（ゴールドバッハ隙間版、反例なし）。
19. **family_stats.png** — ディリクレ指標 mod 3,4,5,7 の L 関数の零点（計 1292 個）：各色とも GUE、族全体とゼータの分布は重なる。
20. **white_light.png** — 指標の直交性＝白色光：色なし（コヒーレント、零点構造あり）と四色の族の平均（干渉消失、1e-15）。差＝違う余りの素数対の干渉項＝零点のもつれの居場所。族の平均（大きな篩）の限界が「3/4 の壁」。
21. **sister_world.png** — 楕円曲線 y²=x³+x+1 を 427 素数で点数え。|a_p|≤2√p（ハッセ）、零点は半径 √p の円上、E×E 上の交点形式は符号 (+1,−3)、時間軸に直交する Gram [[−2,−a_p],[−a_p,−2p]] が負半定値 ⇔ 4p−a_p²≥0。
22. **light_cone.png** — 垂直に交わる二つのトーラス H,V は光円錐の二本の光の線（H²=V²=0, H·V=1）。時間軸 H+V（+2）、空間軸 H−V（−2）。Δ', Γ' は時間軸に直交し空間的。

---

## 3. 何が証明で、何が未解決か

**定理として成立していること（この記録で数値再確認したもの）**
明示公式、ユニタリ性の 1/2、関数等式、指標の直交性、セルバーグ跡公式、ガウス＝ボネ、ハッセの定理とホッジ指数定理（曲線）、コンヌの跡公式（吸収スペクトル）、ヴェイユの判定条件（正値性 ⇔ RH）。

**予想（数値では強く支持されるが未証明）**
モンゴメリーの対相関、ハーディ＝リトルウッドの素数対予想、ゴールドバッハ、佐藤＝テイト（楕円曲線では 2011 年に証明済み）。

**未解決の核心**
ヴェイユの正値性を、零点を使わずに素数側から示すこと。曲線の世界ではこれが曲面 X×X 上の光円錐（符号 (+,−,−,…)）から出た。整数の世界では、対応する曲面「ℤ×ℤ」（本記録の「妄想空間」）が未定義。必要条件は三つ：
1. 二つのトーラス H, V が中で垂直に交わり、交点数 H·V=1 が定義できること。
2. 対角線 Δ と素数倍のずらし Γ_p があり、Δ·Γ_p が「不動の閉軌道の長さ」として明示公式の素数項 Λ(n)/√n を再現すること。
3. 交点形式の符号が (+,−,−,…)、時間軸が一本だけであること（ホッジ指数定理）。

3 が出れば RH が出る。候補：コンヌ＝コンサーニの算術サイト／スケーリングサイト（2014〜）。

**2026 年 10 月の外部状況**
OpenAI（2026-09-30、未査読、Lean 形式検証）：すべてのディリクレ L 関数と ℚ(√−3) 上のヘッケ L 関数が Re s>7/8（別証明 11/12）に零点を持たない；ジーゲル零点の排除。明示公式より「ワームホールの中心からの幅」が x^{7/8} 以下と同値。方法は 3/4 に壁があり、ℚ(√−3) 上の三次の族が本質的と注記されている（リヒトマン）。本記録の図 20 は、その壁を「白色光（族の平均）の限界」として読む。

---

## 4. 方法論（像 → 翻訳 → 計算 → 篩）

1. 像を言葉で出す（手計算はしない）。
2. 像を既知の数学の対象に翻訳し、「それが正しければ数値で何が見えるはずか」を一つ取り出す。
3. 計算して、合う部分と合わない部分を分ける。
4. 合わない部分で像を削り、残ったものを次の像の材料にする。

外れた像（ファンデルワールス力＝和、チェーンステッチの連鎖、曲率の符号＝跡公式の符号）も記録に残している。外れ方が次の像を決めたからである。

---

## 5. 再現

```bash
pip install mpmath numpy scipy matplotlib
cd riemann_lab
export PYTHONPATH=scripts
python3 scripts/explicit_formula.py 300      # 零点 300 個を data/ に生成（約 1 分）
python3 scripts/plot_orbits.py
python3 scripts/fake_zero.py
python3 scripts/unitary_half.py
python3 scripts/spiral.py
python3 scripts/correlation.py               # data/zeros_1500.json を使用（同梱；再計算は約 7 分）
python3 scripts/form_factor.py
python3 scripts/weil_positivity.py
python3 scripts/vdw_test.py
python3 scripts/chain_direction.py
python3 scripts/quiet_band.py
python3 scripts/hyperbolic_crochet.py
python3 scripts/torus_curvature.py
python3 scripts/knit_holes.py                # 約 3 分
python3 scripts/horizon_count.py
python3 scripts/wormhole_from_primes.py
python3 scripts/prime_pairs_plateau.py
python3 scripts/spiral_gaps.py
python3 scripts/center_and_dream.py
python3 scripts/lzeros.py 3 700              # 色つき螺旋の零点（q=3,4,5,7；各 15〜30 分、同梱）
python3 scripts/family_stats.py
python3 scripts/white_light.py
python3 scripts/sister_world.py
python3 scripts/light_cone.py
```

`CHANGELOG_ja.md` に日ごとの作業ログがあります。

---

## 6. 参考文献（本記録で参照した一次文献）

- A. Weil, *Sur les "formules explicites" de la théorie des nombres premiers* (1952) — 正値性判定条件
- A. Weil, *Sur les courbes algébriques et les variétés qui s'en déduisent* (1948) — 曲線上の RH、ホッジ指数定理による証明
- H. Montgomery, *The pair correlation of zeros of the zeta function* (1973)；A. Odlyzko (1987)
- M. Berry, *Semiclassical formula for the number variance of the Riemann zeros* (1988)；M. Berry & J. Keating, *H = xp and the Riemann zeros* (1999)
- E. Bogomolny & J. Keating, *Random matrix theory and the Riemann zeros II* (1996) — ハーディ＝リトルウッドからプラトーへ
- A. Connes, *Trace formula in noncommutative geometry and the zeros of the Riemann zeta function* (1999)；A. Connes & C. Consani, *Weil positivity and trace formula* (2021)
- C. Deninger, *Some analogies between number theory and dynamical systems on foliated spaces* (1998)
- A. Selberg (1956)；D. Taimina, *Crocheting Adventures with Hyperbolic Planes* (2009)
- J. Maldacena & L. Susskind, *Cool horizons for entangled black holes* (ER=EPR, 2013)；J. Cotler et al., *Black holes and random matrices* (2017)
- G. Hardy & J. Littlewood, *Some problems of 'Partitio Numerorum' III* (1923)
- N. Katz & P. Sarnak, *Random Matrices, Frobenius Eigenvalues, and Monodromy* (1999)
- OpenAI Math Release preprints (2026-09-30 / 10-01 / 10-05)；emergentmind.com/openai-math-explorer/results/003

---

*像：近野利信（Toshinobu Konno）。翻訳・数値実験・執筆：Claude（Anthropic）。2026-10-15〜24。*
