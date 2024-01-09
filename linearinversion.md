# Линейная инверсия ( linear inversion)

Прямая задача

## Подготовка данных

#### Модель среды

```mathematica
imp = ConstantArray[1, 51]*2550*2650;
imp[[10 ;; 15]] = 2700 * 2750;
imp[[15 ;; 27]] = 2400 * 2450;
imp[[27 ;; 35]] = 2800 *3000;
ListStepPlot[imp, Filling -> Axis, PlotRange -> All, AspectRatio -> 1/4, ImageSize -> 800, Frame -> True]
```

```mathematica
m = N[(imp[[2 ;;]] - imp[[;; -2]])/(imp[[2 ;;]] + imp[[;; -2]])];
ListPlot[m, Filling -> Axis, PlotRange -> All, AspectRatio -> 1/4, ImageSize -> 800, Frame -> True]
```

#### Модель сигнала

```mathematica
wavelet = Map[WaveletPsi[MexicanHatWavelet[1], #] &, Range[-5, 5, 0.5]];
ListLinePlot[wavelet, Filling -> Axis, PlotRange -> All, AspectRatio -> 1/3, ImageSize -> 500, Frame -> True]
```

## Моделирование сейсмограммы с помощью свёртки

```mathematica
ListConvolve[wavelet, m]
```

```mathematica
ListLinePlot[ListConvolve[wavelet, m], Filling -> Axis, PlotRange -> All, AspectRatio -> 1/4, ImageSize -> 800, Frame -> True]
```

## Моделирование сейсмограммы с помощью умножения матриц

```mathematica
size = Length[m]
```

```mathematica
row = PadRight[wavelet, size]
```

```mathematica
submatrix = Table[row = RotateRight[row], {i, 1, size - Length[wavelet]}];
g = PrependTo[submatrix, PadRight[wavelet, size]];
g // MatrixPlot
```

```mathematica
g // MatrixForm
```

```mathematica
m // MatrixForm
```

```mathematica
d = g . m ;
d // MatrixForm
ListLinePlot[d, Filling -> Axis, PlotRange -> All, AspectRatio -> 1/4, ImageSize -> 800, Frame -> True]
```

Обратная задача

## Инверсия - решение матричного уравнения

#### Матричный метод

![0mzak4zrfe9ob](img\0mzak4zrfe9ob.png)

#### Псевдобратная матрица 

![0dwz81ubg0npi](img\0dwz81ubg0npi.png)

![0s6k3g78yfe4n](img\0s6k3g78yfe4n.png)

#### Вычисление  псевдобратной матрицы

```mathematica
Transpose[g] . Inverse[g . Transpose[g]] // MatrixPlot
```

#### Вычисление столбца коэффициентов отражения (инверсия)

```mathematica
mEstimation = Transpose[g] . Inverse[g . Transpose[g]] . d;
```

```mathematica
Show[ListPlot[m, Filling -> Axis, AspectRatio -> 1/4, ImageSize -> 800, Frame -> True, PlotRange -> All], 
  ListPlot[mEstimation, Filling -> Axis, PlotStyle -> Red, PlotRange -> All]]
```

#### Встроенные способы решения матричных уравнений

```mathematica
mEstimationWolfram = LinearSolve[g, d];
```

## Вычисление акустической жёсткости

![185jdgt9hzjmq](img\185jdgt9hzjmq.png)

```mathematica
impInv = ConstantArray[imp[[1]], Length[imp]];
```

```mathematica
Table[impInv[[i + 1]] = impInv[[i]]*(1 + mEstimation[[i]])/(1 - mEstimation[[i]]), {i, 1, Length[mEstimation] - 1}];
```

```mathematica
Show[ListPlot[imp, Filling -> Axis, PlotRange -> All, AspectRatio -> 1/4, ImageSize -> 800, Frame -> True], 
  ListPlot[impInv, Filling -> Axis, PlotStyle -> Red, PlotRange -> All]]
```

```mathematica
ListLinePlot[d, Filling -> Axis, PlotRange -> All, AspectRatio -> 1/4, ImageSize -> 800, Frame -> True]

```

Выводы

- Устраняется интерференция и появляются границы

- В отличии от условных значений амплитуд отражений появляются  величины напрямую описывающие свойства горных пород (плотность, скорость, акустическая жёсткость)

- Существует неоднозначность в вычислении импеданса - одинаковые коэффициенты отражения могут быть обусловлены разными акустическими свойствами. Это обуславливает необходимость использования скважинной информации