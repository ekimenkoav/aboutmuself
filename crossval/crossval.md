[back to top](./index.html).

# Кроссвалидация при оценке погрешности структурных построений 
[Этот же материал в формате ноутбука Wolfram Mathematica](https://ekimenkoav.github.io/aboutmyself/download/crossvalidation_lecture_notes.nb

***
**Введение**
Рассмотрено использование кроссвалидации в задаче оценки погрешности структурных построений. В "Рекомендациях по использованию материалов сейсморазведки при подсчете запасов" предложено использовать кроссвалидацию, как один из инструментов оценки погрешности прогноза глубин.
Метод заключается в поочередном исключении одной или нескольких скважин из рассмотрения, и использовании их в качестве проверочных.
Наряду с этим подходом традиционно используется и другой способ  - вычисляется стандартное отклонение по набору невязок всех имеющихся скважин. Например с структурная карта построена с использованием скоростей ОГТ или постоянной скоростью.   Тогда  каждая скважина будет характеризоваться невязкой (скорости ОГТ или постоянная скорость не обеспечат точного прогноза глубин). Стандартное отклонение вычисленной по набору невязок и будет мерой погрешности.  
Два описанных способа могут приводить к существенно разным оценкам погрешности.  Это зависит от того как распределены невязки по площади, или другими словами, от того как распределены скоростные аномалии..

**Последовательность вычислений**

1. Подготовка набора ошибок прогноза глубин - *N* точек с координатами *X* и *Y* и величиной ошибки *Z*

1. Выполнение кроссвалидации. 

    - исключение одной точки

    - построение интерполяции по *N-1* точкам

    - вычисление значения *Z* в исключенной точке *delta=Z<sub>true</sub>-Z<sub>interp</sub>*

1. Получение оценки погрешности ( RMS value) по найденному набору значений delta

### Необходимые функции. 2D интерполяция

Ключевым этапом вычислений является построение интерполяции по данным. Здесь используется функция сплайн интерполяции из репозиторий функций Wolfram.

```mathematica
ResourceFunction["PolyharmonicSplineInterpolation"]
```

### Случай некоррелированных ошибок

#### Подготовка данных

Входные данные представляют  50 точек со случайными значениями координат, эти  точки имитируют 50 скважин. В качестве координаты Z будут выступать невязки (error), заданы случайно по нормальному закону со стандартным отклонением заданным переменной sigma.
Фрагмент кода ниже выводит эти точки в виде фрагмента таблицы и в виде гистограммы.

```mathematica
n = 50;
sigma = 20; 
 
xyz = Transpose[{RandomReal[{0, 100}, n], RandomReal[{0, 100}, n], RandomVariate[NormalDistribution[0, sigma], n]}]; 
 
TableForm[Round[xyz[[1 ;; 10]], 0.1], TableHeadings -> {None, {"x", "y", "z"}}]
Histogram[xyz[[All, 3]], {-70, 70, 10}, Frame -> True, ImageSize -> 400, PlotLabel -> "Error distribution", LabelStyle -> Directive[11, Bold, Black]]
```

![фрагмент таблицы](img\cv01_.png)

![гистограмма невязок](img\cv02_.png)



Сначала выполним  интерполяцию, восстановим значения в каждой точке площади и отобразим их на графике.

```mathematica
f = ResourceFunction["PolyharmonicSplineInterpolation"][xyz];
```

```mathematica
Show[
  DiscretePlot3D[f[x, y], {x, 0, 100, 1}, {y, 0, 100, 1}, Filling -> None, Joined -> True, ColorFunction -> Hue], 
  ListPointPlot3D[xyz, PlotRange -> All, PlotStyle -> Black] 
 ]
```

![результат интерполяции](img\cv03_.png)

#### Кроссвалидация

Следующий фрагмент кода выполняет собственно кроссвалидацию.
С помощью функции Table выполняется цикл, в котором итератор  задан от 1 до 50 (Length[xyz]).На каждом шаге исключается проверочная скважина, строится функция интерполяции, восстанавливается значение в проверочной точке. 

Функция StandardDeviation позволяет получить оценки
	как исходного набора невязок StandardDeviation[errors[[All,1]]]
	так и по ошибкам кроссвалидации StandardDeviation[errors[[All,1]]-errors[[All,2]]]

```mathematica
errors = Table[
   	Clear[if]; 
   	tmp = Delete[xyz, i]; 
   	if = ResourceFunction["PolyharmonicSplineInterpolation"][tmp]; 
   	{xyz[[i, 3]], if[xyz[[i, 1]], xyz[[i, 2]]]}, 
    {i, 1, Length[xyz]}];
```

```mathematica
TableForm[Round[{StandardDeviation[errors[[All, 1]]], StandardDeviation[errors[[All, 1]] - errors[[All, 2]]]}, 0.1], TableHeadings -> ]
```

| Метод расчёта | Значение |
| --- | ----------- |
| Simple RMS Value | 20.8 |
| CrossValidation RMS Value | 25.5 |

Выведенная таблица  показывает схожий уровень стандарных отклонений. Т.е. нужно сделать вывод, что в случае невязок заданных случайно оценки по обоим способам дают схожий результат.

### Случай коррелированных ошибок

#### Подготовка данных. Моделирование случайных гауссовых полей.

Удобным способом получения коррелированных по площади данных является моделирование Гауссовых полей. На этой странице используется реализация приведенная по ссылке https://garrettgoon.com/gaussian-fields/.

```mathematica
GaussianRandomField[size : (_Integer?Positive) : 256, dim : (_Integer?Positive) : 2, Pk_ : Function[k, k^-3]] := 
  Module[
    {Pkn, 
    fftIndgen, 
    noise, 
    amplitude, 
    s2}, 
        Pkn = Compile[{{vec, _Real, 1}}, With[{nrm = Norm[vec]}, If[nrm == 0, 0, Sqrt[Pk[nrm]]]], CompilationOptions -> {"InlineExternalDefinitions" -> True}]; 
        s2 = Quotient[size, 2]; 
        fftIndgen = ArrayPad[Range[0, s2], {0, s2 - 1}, "ReflectedNegation"]; 
        noise = Fourier[RandomVariate[NormalDistribution[], ConstantArray[size, dim]]]; 
        amplitude = Outer[Pkn[{##}] &, Sequence @@ ConstantArray[N@fftIndgen, dim]]; 
        InverseFourier[noise*amplitude] 
        ]
```

С помощью моделирования случайного гауссова поля получим карту невязок, в  которой значения коррелируют по площади и меняются плавно. Такая карта может имитировать распределение скоростных аномалий. Например в разрезе присутствует интервал развития соляных куполов, валов штоков. Характер площадного распределения этих объектов и будет обуславливать поведение поля невязок. 

```mathematica
dataRandomField = Re[GaussianRandomField[100]];
MatrixPlot[dataRandomField]
```

![поле случайных коорелированных данных](img\cv04_.png)

#### Подготовка данных. Получение значений в точка скважин.

```mathematica
currentsigma = StandardDeviation[Flatten[dataRandomField]];
coefSigma = sigma/currentsigma;
dataRandomFieldFunction = ListInterpolation[Transpose[Reverse[coefSigma*dataRandomField]]]; (*нужно также указать диапазон интерполяции*)
```

Получим значения невязок в точках скважин

```mathematica
xyzGauss = Table[{x = xyz[[i, 1]], y = xyz[[i, 2]], dataRandomFieldFunction[x, y]}, {i, 1, n}];
Histogram[xyzGauss[[All, 3]], {-70, 70, 10}, Frame -> True, ImageSize -> 400, PlotLabel -> "Error distribution", LabelStyle -> Directive[11, Bold, Black]]
```
![гистограмма случайных коорелированных данных](img\cv05_.png)

Убедимся что смоделированые данные отвечают заранее заданному уровню стандартногот отклонения
```mathematica
StandardDeviation[xyzGauss[[All, 3]]]

(*17.5221*)
```

```mathematica
Clear[ifGauss];
Quiet[ifGauss = ResourceFunction["PolyharmonicSplineInterpolation"][xyzGauss, InterpolationOrder -> 2]];
```

```mathematica

  Show[
   DiscretePlot3D[ifGauss[x, y], {x, 0, 100, 1}, {y, 0, 100, 1}, Filling -> None, Joined -> True, ColorFunction -> Hue], 
   ListPointPlot3D[xyzGauss, PlotRange -> All, PlotStyle -> Black] 
  ]
```
![результат интерполяции сквапжинных данных](img\cv05_.png)


### Кроссвалидация

```mathematica
errorsGauss = Table[
    Clear[ifGauss]; 
    tmp = Delete[xyzGauss, i]; 
    ifGauss = ResourceFunction["PolyharmonicSplineInterpolation"][tmp, InterpolationOrder -> 2]; 
    {xyzGauss[[i, 3]], ifGauss[xyzGauss[[i, 1]], xyzGauss[[i, 2]]]}, 
    {i, 1, Length[xyzGauss]}];
```

```mathematica

  TableForm[Round[{
        StandardDeviation[errorsGauss[[All, 1]]], 
        StandardDeviation[errorsGauss[[All, 1]] - errorsGauss[[All, 2]]]}], 
  ] 
    
   
  
```

| Метод расчёта | Значение |
| --- | ----------- |
| Simple RMS Value | 20 |
| CrossValidation RMS Value | 11 |

---
**Выводы**
* В случае коррелированных по площади ошибок использование кроссвалидации может занизить погрешность. 
* Простой расчёт стандартного отклоения  корректно характеризует погрешность. Но из такой оценки сложно понять тот уровень ошибок кторый можно ожидать в следующей пробуренной скважине. 
* После того как получена струтурная карта, проанализированы ошибки и оценена погрешность карта обязательно приводится в соответсвии с отбивками,а значит в окресности скважин погрешность будет минимальная.В этом случае интересно узнать а на каком расстоянии от скважины погрешность будет достигать своего максимума (значения кторое мы оценили выше). Для этого нужно по скважинным данным оценить степень пространственной корреляции. Об этом подробнее в теме [Вариограмный анализ](https://ekimenkoav.github.io/aboutmyself/index.html)
---

**Ссылки**


1. [Методические рекомендации по использованию данных сейсморазведки (2D, 3D) для подсчета запасов нефти и газа](https://www.geokniga.org/books/12513)

2. [Репозиторий пользовательский функций Wolfram](https://resources.wolframcloud.com/FunctionRepository/)

3. [Детальное описание эффективного алгоритма моделирования случайных гауссовых полей](https://garrettgoon.com/gaussian-fields/)

4. [Обсуждение на stackexchange случайных гауссовых полей](https://mathematica.stackexchange.com/questions/4829/efficiently-generating-n-d-gaussian-random-fields)