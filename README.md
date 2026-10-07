# ideal_flash270

纯 Rust 理想气液平衡闪蒸计算库，无前端、HTTP 服务或联网物性查询。使用 Rust 1.85.1、Serde 1.0.229 和 serde_json 1.0.151。

## 模型与单位

- 温度 `T`：K
- 压力 `P`、饱和蒸气压 `Psat`：bar
- 摩尔分数 `z/x/y`：无因次
- 液体热容 `CpL`、气体热容 `CpV`：J/(mol·K)
- 参考汽化焓 `L` 与焓值：J/mol
- 气相率 `beta`：0 为全液，1 为全气，开区间 `(0,1)` 为两相

Antoine 方程：

```text
log10(Psat) = A - B/(T+C)
```

纯组分理想焓：

```text
hL = CpL*(T-298.15)
hV = L + CpV*(T-298.15)
```

平衡常数为：

```text
K = Psat/P
```

两相 TP 闪蒸使用有界 Rachford-Rice 方程求 `beta`，并计算：

```text
x = z/[1+beta*(K-1)]
y = K*x
h = (1-beta)*sum(x*hL) + beta*sum(y*hV)
```

PH 闪蒸在给定压力和温度括区内对平衡混合焓做二分反求温度。目标必须被端点严格包围；端点恰好命中目标会返回错误，不会把端点当答案。

## 适用范围与校验

- 支持 2 至 16 个组分。
- 组分名必须唯一；所有进料摩尔分数必须为正有限数，且和为 1。
- Antoine 参数、温区、热容和汽化焓必须为有限数。
- `CpL > 0`、`CpV > 0`，压力必须为正有限数。
- 只在所有组分 Antoine 温区的公共交集内计算；公共交集必须有非空内点。
- 公共温区内必须满足 `T+C > 0` 且 `hV-hL > 0`；不做温区外推。
- 全部 `K=1` 的退化状态无法唯一确定相分裂，返回错误。
- 目标未包围、物性越界或迭代未收敛均返回类型化错误。

## 类型化调用

```rust
use ideal_flash270::{Component, Mixture, ph_flash, tp_flash};

let mixture = Mixture::new(vec![
    Component {
        name: "light".into(),
        z: 0.5,
        a: 2.5,
        b: 900.0,
        c: 100.0,
        t_min: 250.0,
        t_max: 450.0,
        cp_l: 100.0,
        cp_v: 80.0,
        latent_heat: 30_000.0,
    },
    // ...其余组分
])?;

let tp = tp_flash(&mixture, 278.0, 1.0)?;
let ph = ph_flash(&mixture, 1.0, 10_000.0, 250.0, 300.0)?;
```

返回结果包含 `temperature`、`pressure`、`phase`、`beta`、按输入顺序排列的 `liquid_composition`/`vapor_composition` 和 `enthalpy`。单相时不存在的相组成返回空数组。

## JSON 调用

`run_json` 接受一个 JSON 对象，用 `kind` 选择 TP 或 PH。参考汽化焓字段名是 `latent_heat`，也接受别名 `L`。

TP 请求：

```json
{
  "kind": "tp",
  "temperature": 278.0,
  "pressure": 1.0,
  "components": [
    {"name":"light","z":0.5,"a":2.5,"b":900.0,"c":100.0,"t_min":250.0,"t_max":450.0,"cp_l":100.0,"cp_v":80.0,"L":30000.0},
    {"name":"heavy","z":0.5,"a":2.5,"b":1000.0,"c":100.0,"t_min":250.0,"t_max":450.0,"cp_l":100.0,"cp_v":80.0,"L":30000.0}
  ]
}
```

PH 请求：

```json
{
  "kind": "ph",
  "pressure": 1.0,
  "enthalpy": 10000.0,
  "t_lo": 250.0,
  "t_hi": 300.0,
  "components": [
    {"name":"light","z":0.5,"A":2.5,"B":900.0,"C":100.0,"t_min":250.0,"t_max":450.0,"CpL":100.0,"CpV":80.0,"L":30000.0},
    {"name":"heavy","z":0.5,"A":2.5,"B":1000.0,"C":100.0,"t_min":250.0,"t_max":450.0,"CpL":100.0,"CpV":80.0,"L":30000.0}
  ]
}
```

字段也接受题面中的大写别名 `T`、`P`、`H`、`A`、`B`、`C`、`Z` 和 `L`，以及 `CpL`、`CpV`。

成功响应中的 `phase` 为 `liquid`、`vapor` 或 `two_phase`。失败时返回 `error` 和机器可读的 `kind`，例如 `invalid_input`、`out_of_range`、`degenerate`、`no_bracket` 或 `not_converged`。

## 构建、测试与示例

```bash
cargo test
cargo run --example tp_ph
cargo build --release
```

项目通过 `.cargo/config.toml` 使用仓库内已准备好的 vendor 依赖，可离线构建。
