# Восстановление значений акустического каротажа по имеющимся кривым ГИС
## 1. Описание задачи
В рамках выполнения работ по интерпретации сейсмических данных решаются задачи прогноза свойств геологического разреза по наблюдениям физических полей.Примеры :
* по измеренным временам прихода нужно прогнозировать глубины,
* по энергии отраженных волн нужно прогнозировать пористость,
* по измеренным значениям акустического каротажа (АК) нужно прогнозировать пористость


Такой же является задача прогноза значений скоростей волн вдоль ствола скважины.Рассмотрим следующую ситуацию.На участке исследования проведены сейсмические работы и пробурено две скважины.В одной из скважин записан полный комплекс ГИС, в т.ч.АК.Во второй скважине запись АК отсутствует.
Для выполнения стратиграфической привязки необходимо во второй скважине восстановить значения АК.

## 2. Загрузка данных по опорной скважине
```
logData=Import[FileNameJoin[{NotebookDirectory[],"dt_predict_11PO_full.xlsx"}]][[1]]
TableForm[logData[[2;;10]],TableHeadings->{None,logData[[1]]}]
```
## 3. Анализ данных
Проанализируем кроссплоты между значениями акустического каротажа и другими каротажными кривыми
```
In[3]:= ListPlot[logData[[2;;,{4,2}]], PlotStyle->Blue, PlotRange->All, Frame->True, ImageSize->500, PlotLabel->"Crossplot NGK vs DT", LabelStyle->Directive[Bold,Orange]]
ListPlot[logData[[2;;,{3,2}]], PlotRange->All, PlotStyle->Blue, Frame->True, ImageSize->500, PlotLabel->"Crossplot GK vs DT", LabelStyle->Directive[Bold,Orange]]
ListPlot[logData[[2;;,{5,2}]], PlotRange->All, PlotStyle->Blue, Frame->True, ImageSize->500, PlotLabel->"Crossplot PZ vs DT", LabelStyle->Directive[Bold,Orange]]
```

## 4. Эмпирические зависимости
В случае когда необходимые каротажные кривые отсутствуют можно воспользоваться эмпирическими зависимостями. Так для восстановления значений АК можно использовать каротаж НГК (формула Заляева) или каротаж сопротивлений (формула Фауста)

````
In[1166]:= TableForm[logData[[2;;10]],TableHeadings->{None,logData[[1]]}]
````
Методика Заляева. Зависимость АК и НГК.
Методика предполагает аппроксимацию наблюденных данных функцией логарифма

```
In[1165]:= ListPlot[logData[[2;;,{4,2}]], PlotStyle->Blue, PlotRange->All, Frame->True, ImageSize->500, PlotLabel->"Crossplot NGK vs DT", LabelStyle->Directive[Bold,Orange]]
DT =a*Ln[NGK+b]+c, 
```

DT 		время пробега 
NGK	 	значения кривой нейтронного каротажа
a, b, c 	параметры модели, которые нужно оценить
Функция FindFit позволяет задать модель для аппроксимации данных (формулу) и оценить коэффициенты этой модели. В данном примере используется натуральный логарифм.
``
		FindFit[data, a*Log[x + b] + c, {a, b, c}, {x}]
In[1179]:= Clear[a,b,c]
FindFit[logData[[2;;,{4,2}]],a*Log[x+b]+c,{a,b,c},{x}]
``
Получив искомые значения параметров a, b и c можно создать функцию для прогноза значений акустического каротажа
``
In[1181]:= fNGK[x_]:=-27.797*Log[x-1.57585]+202.99
In[1183]:= 
Show[
ListPlot[logData[[2;;;;10,{4,2}]],PlotRange->{{0,10}, {120,350}},PlotStyle->Blue, Frame->True],
Plot[fNGK[x],{x,0,10},PlotStyle->Red,PlotRange->{{0,10}, {120,350}}]
]
dtNGK=fNGK[logData[[2;;, 4]]]
In[1192]:= Show[ListLinePlot[logData[[2;;,{1,2}]], PlotStyle->Black, PlotRange->{All,{120,350}}, AspectRatio->1/6, ImageSize->700],ListLinePlot[Transpose[{logData[[2;;,1]], dtNGK}], PlotStyle->Blue, PlotRange->{All,{120,350}}], PlotLabel->"DT by NGK (Zalyaev)", LabelStyle->Directive[Bold,Orange]]
``

Формула Фауста. Зависимость АК и каротажа сопротивлений
DT =a*(R)^b, 

DT - 		время пробега 
R - 		сопротивления
a, b -		параметры модели, которые нужно оценить