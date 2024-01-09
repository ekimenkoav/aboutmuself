# Восстановление значений акустического каротажа по имеющимся кривым ГИС

## 1. Описание задачи

В рамках выполнения работ по интерпретации сейсмических данных решаются задачи прогноза свойств геологического разреза по наблюдениям физических полей.Примеры :
*по измеренным временам прихода нужно прогнозировать глубины,
*по энергии отраженных волн нужно прогнозировать пористость,
*по измеренным значениям акустического каротажа (АК) нужно прогнозировать пористость

Такой же является задача прогноза значений скоростей волн вдоль ствола скважины.Рассмотрим следующую ситуацию.На участке исследования проведены сейсмические работы и пробурено две скважины.В одной из скважин записан полный комплекс ГИС, в т.ч.АК.Во второй скважине запись АК отсутствует.
Для выполнения стратиграфической привязки необходимо во второй скважине восстановить значения АК.

## 2. Загрузка данных по опорной скважине

```mathematica
logData = Import[FileNameJoin[{NotebookDirectory[], "dt_predict_11PO_full.xlsx"}]][[1]]
TableForm[logData[[2 ;; 10]], TableHeadings -> {None, logData[[1]]}]
```

## 3. Анализ данных

Проанализируем кроссплоты между значениями акустического каротажа и другими каротажными кривыми

```mathematica
ListPlot[logData[[2 ;;, {4, 2}]], PlotStyle -> Blue, PlotRange -> All, Frame -> True, ImageSize -> 500, PlotLabel -> "Crossplot NGK vs DT", LabelStyle -> Directive[Bold, Orange]]
ListPlot[logData[[2 ;;, {3, 2}]], PlotRange -> All, PlotStyle -> Blue, Frame -> True, ImageSize -> 500, PlotLabel -> "Crossplot GK vs DT", LabelStyle -> Directive[Bold, Orange]]
ListPlot[logData[[2 ;;, {5, 2}]], PlotRange -> All, PlotStyle -> Blue, Frame -> True, ImageSize -> 500, PlotLabel -> "Crossplot PZ vs DT", LabelStyle -> Directive[Bold, Orange]]
```

## 4. Эмпирические зависимости

В случае когда необходимые каротажные кривые отсутствуют можно воспользоваться эмпирическими зависимостями. Так для восстановления значений АК можно использовать каротаж НГК (формула Заляева) или каротаж сопротивлений (формула Фауста)

```mathematica
TableForm[logData[[2 ;; 10]], TableHeadings -> {None, logData[[1]]}]
```

### Методика Заляева. Зависимость АК и НГК.

Методика предполагает аппроксимацию наблюденных данных функцией логарифма

```mathematica
ListPlot[logData[[2 ;;, {4, 2}]], PlotStyle -> Blue, PlotRange -> All, Frame -> True, ImageSize -> 500, PlotLabel -> "Crossplot NGK vs DT", LabelStyle -> Directive[Bold, Orange]]
```

DT =a*Ln[NGK+b]+c, 

DT 		время пробега 
NGK	 	значения кривой нейтронного каротажа
a, b, c 	параметры модели, которые нужно оценить

Функция FindFit позволяет задать модель для аппроксимации данных (формулу) и оценить коэффициенты этой модели. В данном примере используется натуральный логарифм.

		FindFit[data, a*Log[x + b] + c, {a, b, c}, {x}]

```mathematica
Clear[a, b, c]
FindFit[logData[[2 ;;, {4, 2}]], a*Log[x + b] + c, {a, b, c}, {x}]
```

Получив искомые значения параметров a, b и c можно создать функцию для прогноза значений акустического каротажа

```mathematica
fNGK[x_] := -27.797*Log[x - 1.57585] + 202.99
```

```mathematica

  Show[
   ListPlot[logData[[2 ;; ;; 10, {4, 2}]], PlotRange -> { {0, 10}, {120, 350} }, PlotStyle -> Blue, Frame -> True], 
   Plot[fNGK[x], {x, 0, 10}, PlotStyle -> Red, PlotRange -> {{0, 10}, {120, 350}}] 
  ]
```

```mathematica
dtNGK = fNGK[logData[[2 ;;, 4]]]
```

```mathematica
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtNGK}], PlotStyle -> Blue, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by NGK (Zalyaev)", LabelStyle -> Directive[Bold, Orange]]
```

### Формула Фауста. Зависимость АК и каротажа сопротивлений

DT =a*$(R)^b$, 

DT - 		время пробега 
R - 		сопротивления
a, b -		параметры модели, которые нужно оценить

```mathematica
Clear[a, b]
FindFit[logData[[2 ;;, {1, 5, 2}]], a*(x*y)^b, {a, b}, {x, y}]
```

```mathematica
(*fR[x_]:=234.748*(x)^(-0.0690263)*)
  fR[{x_, y_}] := 319.964*(x*y)^(-0.0494973)
```

```mathematica
Show[
  ListPlot[logData[[2 ;; ;; 10, {5, 2}]], PlotRange -> {{0, 250}, {120, 350}}, PlotStyle -> Blue, Frame -> True], 
  Plot[fR[x], {x, 0, 250}, PlotStyle -> Red, PlotRange -> {{0, 250}, {120, 350}}] 
 ]
```

```mathematica
dtFaust = Map[fR, logData[[2 ;;, {1, 5}]]]
```

```mathematica
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtFaust}], PlotStyle -> Red, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by Faust", LabelStyle -> Directive[Bold, Orange]]
```

## Методы машинного обучения

Язык Wolfram Language включает в себя широкий спектр современных  возможностей машинного обучения таких как прогнозирование и классификация.
Функции работают со многими типами данных, включая числовые, категориальные, временные ряды, текстовые, изображения и аудио.

В данном примере в качестве данных данных выступают каротажи. Задача заключается в том, чтобы по измеренным замерам каротажей  "GR",  "NGK" и  "PZ" прогнозировать "DT"

```mathematica
TableForm[logData[[2 ;; 10]], TableHeadings -> {None, logData[[1]]}]
```

Встроенная функция **Predict** является основным инструментом для прогнозирования с использованием всех методов машинного обучения: 

**"DecisionTree"** 			predict using a decision tree
**"GradientBoostedTrees"**	predict using an ensemble of trees trained with gradient boosting
**"LinearRegression"**		predict from linear combinations of features
**"NearestNeighbors"**		predict from nearest neighboring examples
**"NeuralNetwork"**			predict using an artificial neural network
"**RandomForest"**			predict from Breiman\[Dash]Cutler ensembles of decision trees
**"GaussianProcess"**		predict using a Gaussian process prior over functions

### Обучающая выборка

```mathematica
forPredict = Map[#[[{1, 4, 5}]] -> #[[2]] &, logData[[2 ;;]]];
TableForm[forPredict[[1 ;; 10]]]
```

### Множественная регрессия ([ссылка](https://www.youtube.com/watch?v=zzoVdPgVofM&ysclid=lp80o92v6b922626680))

```mathematica
predictLM = Predict[forPredict, Method -> "LinearRegression"]
```

```mathematica
PredictorInformation[predictLM, "Function"]
```

```mathematica
dtLM = predictLM[forPredict[[All, 1]]]
```

```mathematica
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtLM}], PlotStyle -> Green, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by Linear regression model", LabelStyle -> Directive[Bold, Orange]]
```

### Случайный лес ([ссылка](https://www.youtube.com/watch?v=nbxiRdAk1JY&ysclid=lp80lj5aub913635003))

```mathematica
predictRF = Predict[forPredict, Method -> "RandomForest"]
```

```mathematica
dtRF = predictRF[forPredict[[All, 1]]]
```

```mathematica
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtRF}], PlotStyle -> Orange, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by Random forest", LabelStyle -> Directive[Bold, Orange]]
```

### Нейронные сети

```mathematica
predictNN = Predict[forPredict, Method -> "NeuralNetwork"]
```

```mathematica
dtNN = predictNN[forPredict[[All, 1]]]
```

```mathematica
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtNN}], PlotStyle -> Magenta, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by Neural Network", LabelStyle -> Directive[Bold, Orange]]
```

## Сопоставление результатов

```mathematica
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtFaust}], PlotStyle -> Red, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by Faust", LabelStyle -> Directive[Bold, Orange]]
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtNGK}], PlotStyle -> Blue, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by NGK (Zalyaev)", LabelStyle -> Directive[Bold, Orange]]
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtLM}], PlotStyle -> Green, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by Linear regression model", LabelStyle -> Directive[Bold, Orange]]
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtRF}], PlotStyle -> Orange, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by Random forest", LabelStyle -> Directive[Bold, Orange]]
Show[ListLinePlot[logData[[2 ;;, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logData[[2 ;;, 1]], dtNN}], PlotStyle -> Magenta, PlotRange -> {All, {120, 350}}], PlotLabel -> "DT by Neural Network", LabelStyle -> Directive[Bold, Orange]]
```

## Использование моделей прогноза для новой скважины

### Каротажи по новой скважине

```mathematica
logDataTest = Import[FileNameJoin[{NotebookDirectory[], "dt_predict_2PO_noDT.xlsx"}]][[1]];
TableForm[logDataTest[[2 ;; 10]], TableHeadings -> {None, logDataTest[[1]]}]
```

### Расчёт кривой АК с использованием ранее подготовленных моделей 

```mathematica
dtTestNN = predictNN[logDataTest[[2 ;;, {1, 3, 4}]]]
dtTestRF = predictRF[logDataTest[[2 ;;, {1, 3, 4}]]]
dtTestLM = predictLM[logDataTest[[2 ;;, {1, 3, 4}]]]
dtTestNGK = fNGK[logDataTest[[2 ;;, 3]]]
dtTestFaust = fR[logDataTest[[2 ;;, 4]]]
```

### Контроль

```mathematica
logdataControl = Import[FileNameJoin[{NotebookDirectory[], "dt_predict_2PO_full.xlsx"}]][[1]];
TableForm[logdataControl[[2 ;; 10, {1, 2}]], TableHeadings -> {None, logdataControl[[1, {1, 2}]]}]
```

```mathematica
errorRmsNGK = Round[StandardDeviation[logdataControl[[2 ;; -500, 2]] - dtTestNGK[[1 ;; -500]]]];
errorRmsLM = Round[StandardDeviation[logdataControl[[2 ;; -500, 2]] - dtTestLM[[1 ;; -500]]]];
errorRmsRF = Round[StandardDeviation[logdataControl[[2 ;; -500, 2]] - dtTestRF[[1 ;; -500]]]];
errorRmsNN = Round[StandardDeviation[logdataControl[[2 ;; -500, 2]] - dtTestNN[[1 ;; -500]]]];
errorRmsFaust = Round[StandardDeviation[logdataControl[[2 ;; -500, 2]] - dtTestFaust[[1 ;; -500]]]];
```

```mathematica
Show[ListLinePlot[logdataControl[[2 ;; -500, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logDataTest[[2 ;; -500, 1]], dtTestFaust[[1 ;; -500]]}], PlotStyle -> Red, PlotRange -> {All, {120, 350}}], PlotLabel -> StringJoin["DT by Faust, RMS error is ", ToString[errorRmsFaust], " ms"], LabelStyle -> Directive[Bold, Orange]]
Show[ListLinePlot[logdataControl[[2 ;; -500, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logDataTest[[2 ;; -500, 1]], dtTestNGK[[1 ;; -500]]}], PlotStyle -> Blue, PlotRange -> {All, {120, 350}}], PlotLabel -> StringJoin["DT by NGK (Zalyaev), RMS error is ", ToString[errorRmsNGK], " ms"], LabelStyle -> Directive[Bold, Orange]]
Show[ListLinePlot[logdataControl[[2 ;; -500, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logDataTest[[2 ;; -500, 1]], dtTestLM[[1 ;; -500]]}], PlotStyle -> Green, PlotRange -> {All, {120, 350}}], PlotLabel -> StringJoin["DT by Linear regression model is ", ToString[errorRmsLM], " ms"], LabelStyle -> Directive[Bold, Orange]]
Show[ListLinePlot[logdataControl[[2 ;; -500, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logDataTest[[2 ;; -500, 1]], dtTestRF[[1 ;; -500]]}], PlotStyle -> Orange, PlotRange -> {All, {120, 350}}], PlotLabel -> StringJoin["DT by Random forest is ", ToString[errorRmsRF], " ms"], LabelStyle -> Directive[Bold, Orange]]
Show[ListLinePlot[logdataControl[[2 ;; -500, {1, 2}]], PlotStyle -> Black, PlotRange -> {All, {120, 350}}, AspectRatio -> 1/6, ImageSize -> 700], ListLinePlot[Transpose[{logDataTest[[2 ;; -500, 1]], dtTestNN[[1 ;; -500]]}], PlotStyle -> Magenta, PlotRange -> {All, {120, 350}}], PlotLabel -> StringJoin["DT by Neural Network is ", ToString[errorRmsNN], " ms"], LabelStyle -> Directive[Bold, Orange]]
```

