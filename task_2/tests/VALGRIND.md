# Valgrind на CachyOS

Установленные сборки glibc и Valgrind для x86-64-v4 содержат AVX-512.
Valgrind 3.25.1 падает с SIGILL ещё в загрузчике или собственной библиотеке
`vgpreload_memcheck`. Системные стартовые объектные файлы также добавляют
исполняемому файлу требование x86-64-v4, даже при обычном `-march=x86-64`.

Запуск с совместимыми библиотеками и стартовыми объектными файлами:

```sh
make valgrind
make valgrind VALGRIND_ARGS='5 2 5 1'
make valgrind VALGRIND_ARGS='3 2 3 0 test_data/identity_3.txt'
```

Первый запуск загружает стандартные пакеты glibc, libstdc++, libgcc и Valgrind
из Arch Linux в `build/`. Нужны `curl`, `bsdtar` и доступ к интернету.
Для отладочных символов загрузчика нужен `debuginfod-find` (пакет debuginfod)
и переменная окружения:

```sh
export DEBUGINFOD_URLS=https://debuginfod.archlinux.org
```

В CachyOS она обычно уже настроена. Valgrind загружает символы автоматически.
Если появляется `Fatal error at startup` с упоминанием `memcmp`, проверьте
доступность debuginfod и эту переменную.

Сборка `build/a.valgrind` использует те же объектные файлы проекта, но
совместимые стартовые объекты glibc. Системные пакеты и обычный `a.out`
не заменяются. Ошибки памяти и утечки дают код возврата 99.
После `make clean` окружение будет загружено заново.
