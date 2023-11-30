# Bachelor-Thesis
Код программы для среды MATLAB, легшей в основу моей выпускной квалификационной работы по специальности "Оптотехника".

Тема: "Гиперспектрометр на основе мультиспектральной камеры"

Целью ставилось проведение математического моделирования, воспроизводившего процесс последовательной съемки картинки обычной трехканальной RGB-камерой через набор оптических светофильтров.
Благодаря полученному массиву сигналов напряжения, методами решения обратной задачи полностью воссоздавался оригинальный сигнал спектральной плотности яркости изборажения.
Такой подход обосновывался тем, что приборы с высокой спектральной разрешающей способностью - гиперспектрометры, позволяющие получить информативный снимок спектра объекта - во-первых, достаточно сложны в производстве и дорогостоящи, а во-вторых, обладают более низким пространственным разрешением по сравнению с RGB-камерами. 

#### ***DISCLAIMER*** ####
***Алгоритм был написан в далёком 2019 году, писался с нуля в незнакомой на тот момент среде, без какой-либо оптимизации кода, без ориентировки на читабельность и легкость внесения изменений (многие константы не были прописаны и задавались вручную, некоторые повторяющиеся операции скорее всего можно было свести к более компактным функциям).***

В репозитории прилагаются файлы: 
1) "algorithm.m", содержащий основную программу
2) "pic2val.m", содержащий функцию для перевода "сигнала яркости"/"кривой пропускания светофильтров" в массив точек для последующей обработки
3) .png картинки с вышеупомянутыми кривыми яркости объектов и коэффициентов пропускания выбранных светофильтров, используемых для моделирования процесса съемки
