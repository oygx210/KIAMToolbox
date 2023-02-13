# KIAM Astrodynamics Toolbox

KIAM Astrodynamics Toolbox - это библиотека, написанная на языках Fortran и Python в Институте прикладной математики им. М.В. Келдыша.
Библиотека содержит астрономические постоянные, параметры небесных тел, правые части уравнений движения в различных моделях, методы интегрирования обыкновенных дифференциальных уравнений, функции преобразования переменных, систем координат и единиц измерения, а также высокоуровневый класс для проектирования траекторий.

## О версиях

Текущая версия: 1.4.

Следующая версия 1.5 выйдет 13.02.2023.

## Содержимое библиотеки

Спустя какое-то время здесь появится справка по использованию библиотеки. Пока что примеры можно посмотреть в файле kiam_examples.py.
Коротко опишу содержимое файлов:

- `FKIAMToolbox.cp39-win_amd64` – файл со скомпилированными с использованием F2PY на языке Fortran астродинамическими функциями. Предполагается, что пользователи не будут его изменять и импортировать в своих файлах.

- `kiam.py` – основной файл, содержащий функции-интерфейсы для обращения к файлу выше, именно с ним предполагается основная работа пользователей. Таким образом, обращение к низкоуровневым функциям осуществляется через kiam.py, а не напрямую к FKIAMToolbox.cp39-win_amd64.

- `kiam_examples.py` – файл, который содержит несколько примеров использования библиотеки.

- `Trajectory.py` – высокоуровневый класс для проектирования траекторий, он облегчает распространение траекторий, перевод между системами координат и замены переменных.

- `JPLEPH` – файл с эфемеридами небесных тел, используется функциями из FKIAMToolbox.cp39-win_amd64, разработанные в NASA JPL. Это файл содержит эфемериды D430, которые можно также скачать самостоятельно с сайта https://ssd.jpl.nasa.gov/.

- `dll-файлы` – вспомогательные файлы для работы Fortran-Python-интерфейса из Intel(r) Visual Fortran Compiler, которые также можно скачать отдельно.

## Системные требования

ОС Windows 10.

См. также файл `requirements.txt`, в котором перечислены пакеты-зависимости и их версии.

## Информация и условия использования

Справка по модулю kiam.py: https://shmaxg.github.io/KIAMToolbox/html/kiam.html

Справка по модулю Trajectory: https://shmaxg.github.io/KIAMToolbox/html/Trajectory.html

Библиотека активно разрабатывается и улучшается, в ближайший год-два возможно существенное изменение каких-то ее частей.

Библиотека свободна для использования (используется MIT License), единственное просьба - упоминать ее использование в благодарностях в научных трудах, например так:

"Библиотека KIAM Astrodynamics Toolbox, Широбоков М.Г., MIT license, см. https://github.com/shmaxg/KIAMToolbox"

"Package KIAM Astrodynamics Toolbox by Maksim Shirobokov, MIT license, see https://github.com/shmaxg/KIAMToolbox for details"

В дальнейшем здесь появится ссылка на статью, описывающая библиотеку, на которую также можно ссылаться.

По различным вопросам можно обращаться к основному автору библиотеки Широбокову Максиму Геннадьевичу: shirobokov@keldysh.ru.
