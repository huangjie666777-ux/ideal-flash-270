# ideal_flash270

理想体系气液闪蒸计算库（Rust 1.85.1，serde 1.0.229 / serde_json 1.0.151，无网络依赖）。
用于流程设计中判断混合进料经减压或换热后的气液分配：支持 **TP 闪蒸**（给定温压求平衡）
与 **PH 闪蒸**（给定压力与目标摩尔焓反求平衡温度）。

## 单位与模型

- 温度 K，压力 bar，焓 J/mol。
- 饱和蒸气压：`log10(Psat/bar) = A - B / (T + C)`，仅在各自有效温区 `[t_min, t_max]` 内使用，**不外推**。
- 液相摩尔焓：`hL = CpL * (T - 298.15)`；气相摩尔焓：`hV = L + CpV * (T - 298.15)`。
- `CpL`、`CpV` 为正的常数比热，`L` 为 298.15 K 参考汽化焓（正）。

## 适用范围与校验

- 2 至 16 个组分，名称唯一；所有 `z > 0` 且总和为 1；拒绝任何非有限数值（NaN/Inf）。
- 公共温区为各组分 Antoine 有效温区的交集；在公共温区内要求 `T + C > 0` 且 `hV - hL > 0`，否则拒绝建模。
- TP：K 因子 `K_i = Psat_i / P`，用 Rachford-Rice 判相并在 `[0,1]` 上有界二分求气相率 β。
  两相时 `x_i = z_i / (1 + β(K_i - 1))`、`y_i = K_i x_i`，保证组分守恒与各相归一。
  单相时 β 为 0（液）或 1（气），存在相组成为 z，缺失相返回空数组。
  全部 `K_i = 1` 的退化状态明确报错（`DegenerateEquilibrium`）。
- PH：以平衡态混合焓 `h(T)` 对温度二分求解，允许跨越单相与两相区；
  目标焓未被严格包围、物性越界或迭代不收敛均返回错误，**端点不作为答案**。

## 代码结构

- `src/props.rs` — 组分数据、物性校验、Psat/焓计算
- `src/equilibrium.rs` — Rachford-Rice 相平衡（TP 闪蒸）
- `src/energy.rs` — PH 能量求解
- `src/serde_io.rs` — JSON 输入输出
- `tests/flash_tests.rs` — 自测；`examples/flash_demo.rs` — 可运行示例

## 调用方式

类型化 API：

```rust
use ideal_flash270::props::{validate_system, Component};
use ideal_flash270::{tp_flash, ph_flash};

let sys = validate_system(components)?;          // Vec<Component>
let r = tp_flash(&sys, 370.0, 1.01325)?;         // T [K], P [bar]
let r2 = ph_flash(&sys, 1.01325, 17000.0, 280.0, 410.0)?; // P, h [J/mol], 温度括区
```

JSON API（`flash_tp_json` / `flash_ph_json`，输入输出均为 JSON 字符串）：

```json
{
  "components": [
    {"name": "benzene", "z": 0.4,
     "antoine": {"a": 4.01814, "b": 1203.835, "c": -53.226, "t_min": 280.0, "t_max": 380.0},
     "cp_l": 136.0, "cp_v": 82.0, "l_ref": 30760.0},
    {"name": "toluene", "z": 0.6,
     "antoine": {"a": 4.07827, "b": 1343.943, "c": -53.773, "t_min": 280.0, "t_max": 410.0},
     "cp_l": 157.0, "cp_v": 103.0, "l_ref": 33180.0}
  ],
  "temperature": 370.0,
  "pressure": 1.01325
}
```

PH 请求将 `temperature` 换成 `enthalpy`、`t_lo`、`t_hi`。结果包含
`temperature`、`pressure`、`phase`（`liquid` / `vapor` / `two_phase`）、`beta`、
按输入顺序排列的 `x` / `y` 及混合焓 `enthalpy`；失败时返回 `{"error": "..."}`。

## 构建与运行

```sh
cargo test --offline              # 自测
cargo run --offline --example flash_demo   # 苯/甲苯 TP + PH 示例
```
