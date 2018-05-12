# Aerodnamics Characteristics Calculation of MultiRotors Using Various Models
---
## Rotor Parameters Moudle

## Uniform Inflow Moudle

## Rotor Force Moudle

## Flapping Moudle
### 隐式迭代求解步骤
1. 初始化旋翼参数
2. 初始化均匀入流场，给定挥舞状态量初值
3. 根据入流和挥舞状态量求解叶素气动力，根据叶素气动力积分求解旋翼力和力矩
4. 判断旋翼力和力矩是否平衡，若平衡，则输出
5. 若不平衡，根据旋翼力和力矩，代入挥舞函数，反向求解挥舞状态量，并更新操纵输入量，
6. 将新的挥舞状态量带回步骤3重新计算
