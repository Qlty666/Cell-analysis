# receptors

把受体文件放到这里，默认名称为 `receptor.pdb`。

- 支持 PDB / mol2 / PDBQT 输入。
- 如果是 PDBQT，系统会直接复制为 `receptor.pdbqt`。
- 对接盒中心和尺寸在 `config/docking_config.json` 中设置。
- 建议先移除水分子和共结晶配体，只保留蛋白链。
