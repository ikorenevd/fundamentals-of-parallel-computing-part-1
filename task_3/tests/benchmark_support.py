#!/usr/bin/env python3
"""Запуск без ограничения времени, чтение отчёта и добавление строки CSV."""

import argparse
import csv
import fcntl
import io
import json
import math
import os
from pathlib import Path
import platform
import re
import shlex
import signal
import subprocess
import sys


REPORT = re.compile(
    r".+ : Task = (?P<task>\d+) Res1 = (?P<res1>\S+)"
    r" Res2 = (?P<res2>\S+) T1 = (?P<t1>\S+) T2 = (?P<t2>\S+)"
    r" S = (?P<s>-?\d+) N = (?P<n>-?\d+) M = (?P<m>-?\d+)"
)
EXPECTATIONS = {"solved", "unsolved", "error", "any"}
SYSTEM_COLUMNS = ("system", "architecture", "cpu", "hostname")


def collect_system_info():
    """Сведения о машине запуска; недоступные поля обозначаются unknown."""
    system = platform.platform(aliased=True)
    cpu = platform.processor()
    if platform.system() == "Linux":
        try:
            release = {}
            for line in Path("/etc/os-release").read_text(encoding="utf-8").splitlines():
                key, separator, value = line.partition("=")
                if separator:
                    words = shlex.split(value)
                    release[key] = words[0] if words else ""
            if release.get("PRETTY_NAME"):
                system = f"{release['PRETTY_NAME']} (Linux {platform.release()})"
        except (OSError, ValueError):
            pass
        try:
            entries = {}
            for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
                key, separator, value = line.partition(":")
                if separator and value.strip():
                    entries.setdefault(key.strip(), value.strip())
            cpu = entries.get("model name") or entries.get("Hardware") or cpu
        except OSError:
            pass
    elif platform.system() == "Darwin":
        try:
            result = subprocess.run(["/usr/sbin/sysctl", "-n", "machdep.cpu.brand_string"],
                                    capture_output=True, text=True, check=False)
            if result.returncode == 0 and result.stdout.strip():
                cpu = result.stdout.strip()
            else:
                cpu = ""
        except OSError:
            cpu = ""
    info = dict(system=system, architecture=platform.machine(),
                cpu=cpu, hostname=platform.node())
    return {key: " ".join(value.split()) or "unknown" for key, value in info.items()}


def _validate_system_info(info):
    if not isinstance(info, dict) or set(info) != set(SYSTEM_COLUMNS):
        raise ValueError("сведения о системе должны содержать system, architecture, cpu, hostname")
    for key in SYSTEM_COLUMNS:
        value = info[key]
        if not isinstance(value, str) or not value.strip() or "\n" in value or "\r" in value:
            raise ValueError(f"некорректные сведения о системе: {key}")
    return info


def _parse_summary(stdout, n, m, s):
    lines = stdout.rstrip().splitlines()
    report = REPORT.fullmatch(lines[-1]) if lines else None
    if report is None:
        raise ValueError("нет корректной итоговой строки")
    values = report.groupdict()
    if tuple(int(values[key]) for key in ("task", "n", "m", "s")) != (11, n, m, s):
        raise ValueError("параметры итоговой строки не совпадают с тестом")
    numbers = {key: float(values[key]) for key in ("res1", "res2", "t1", "t2")}
    if not all(math.isfinite(value) for value in numbers.values()):
        raise ValueError("итоговая строка содержит NaN или бесконечность")
    if numbers["t1"] < 0 or numbers["t2"] < 0:
        raise ValueError("отрицательное время в итоговой строке")
    if numbers["res1"] == -1 and numbers["res2"] == -1:
        status = "unsolved"
    elif numbers["res1"] >= 0 and numbers["res2"] >= 0:
        status = "solved"
    else:
        raise ValueError("некорректные отрицательные невязки")
    return values, numbers, status


def parse_report(stdout, n, m, s):
    values, _numbers, status = _parse_summary(stdout, n, m, s)
    if status == "unsolved":
        return "-1", "-1", status
    return values["t1"], values["t2"], status


def _kill_group(process):
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def execute_case(command, log_path, expectation, n, m, s):
    if expectation not in EXPECTATIONS:
        raise ValueError("неизвестное ожидание теста")
    stdout, stderr, detail = "", "", ""
    t1, t2, status, passed, returncode = "-1", "-1", "error", False, None
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True,
                                   encoding="utf-8", errors="replace",
                                   start_new_session=True)
    except OSError as error:
        detail = str(error)
    else:
        try:
            stdout, stderr = process.communicate()
        except BaseException:
            _kill_group(process)
            process.communicate()
            raise
        else:
            if process.returncode != 0:
                status = "crash" if process.returncode < 0 else "error"
                passed = status == "error" and expectation == "error"
                detail = f"код завершения {process.returncode}"
            else:
                try:
                    values, numbers, status = _parse_summary(stdout, n, m, s)
                except ValueError as error:
                    status, detail = "invalid_report", str(error)
                else:
                    if status == "solved":
                        t1, t2 = values["t1"], values["t2"]
                    passed = expectation == status or expectation == "any"
                    # Для матриц Гильберта значения невязок диагностические:
                    # общий допуск не учитывает их плохую обусловленность.
                    if expectation == "solved" and status == "solved" and s != 4:
                        passed = numbers["res1"] <= 1e-8 and numbers["res2"] <= 1e-6
                        if not passed:
                            detail = "невязки превышают Res1=1e-8 или Res2=1e-6"
                    if not passed and not detail:
                        detail = f"ожидалось {expectation}, получено {status}"
        returncode = process.returncode
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        f"Command: {shlex.join(map(str, command))}\nExit: {returncode}\n"
        f"Status: {status}\nExpectation: {expectation}\nPassed: {passed}\n"
        f"Detail: {detail}\n\n--- stdout ---\n{stdout}"
        f"\n--- stderr ---\n{stderr}", encoding="utf-8")
    return t1, t2, status, passed


def read_manifest(path):
    cases, names = [], set()
    with Path(path).open(encoding="utf-8", newline="") as source:
        for line, fields in enumerate(csv.reader(source, delimiter="\t"), 1):
            if len(fields) != 7:
                raise ValueError(f"manifest: строка {line}: ожидалось 7 полей")
            name, n, m, r, s, filename, expectation = fields
            if not re.fullmatch(r"[A-Za-z0-9_-]+", name) or name in names:
                raise ValueError(f"manifest: некорректное или повторное имя {name!r}")
            if expectation not in EXPECTATIONS or not filename:
                raise ValueError(f"manifest: некорректный тест {name}")
            cases.append(dict(name=name, n=int(n), m=int(m), r=int(r), s=int(s),
                              filename=filename, expectation=expectation))
            names.add(name)
    if not cases:
        raise ValueError("manifest: нет тестов")
    return cases


def result_header(cases):
    return ["timestamp", *SYSTEM_COLUMNS] + [f"{case['name']}.{metric}"
                                             for case in cases for metric in ("T1", "T2")]


def _read_open_results(target, header):
    target.seek(0)
    contents = target.read()
    if not contents:
        return False, False, [], contents
    if not contents.endswith("\n"):
        raise ValueError("файл результатов содержит незавершённую строку")
    rows = csv.reader(io.StringIO(contents), strict=True)
    stored_header = next(rows, None)
    legacy_header = ["timestamp", *header[1 + len(SYSTEM_COLUMNS):]]
    legacy = stored_header == legacy_header
    if stored_header != header and not legacy:
        raise ValueError("набор тестов отличается от заголовка CSV; выберите другой --results")
    offset = 1 if legacy else 1 + len(SYSTEM_COLUMNS)
    history = []
    for line, row in enumerate(rows, 2):
        if len(row) != len(stored_header) or not row[0].strip():
            raise ValueError(f"CSV: некорректная строка {line}")
        if not legacy:
            _validate_system_info(dict(zip(SYSTEM_COLUMNS, row[1:offset])))
        numbers = []
        for value in row[offset:]:
            number = float(value)
            if not math.isfinite(number) or (number < 0 and number != -1):
                raise ValueError(f"CSV: некорректное время в строке {line}")
            numbers.append(number)
        for index in range(0, len(numbers), 2):
            if (numbers[index] == -1) != (numbers[index + 1] == -1):
                raise ValueError(f"CSV: неполная пара отказа T1/T2 в строке {line}")
        history.append(row)
    return True, legacy, history, contents


def check_results(results, cases):
    results = Path(results)
    if results.exists():
        with results.open(encoding="utf-8", newline="") as target:
            fcntl.flock(target, fcntl.LOCK_SH)
            _read_open_results(target, result_header(cases))


def append_results(results, cases, measurements, timestamp, system_info=None):
    if not timestamp.strip() or "\n" in timestamp or "\r" in timestamp:
        raise ValueError("некорректная дата запуска")
    with Path(measurements).open(encoding="utf-8", newline="") as source:
        rows = list(csv.reader(source, delimiter="\t"))
    if len(rows) != len(cases):
        raise ValueError("число результатов не совпадает с числом тестов")
    info = _validate_system_info(collect_system_info() if system_info is None else system_info)
    row = [timestamp, *(info[key] for key in SYSTEM_COLUMNS)]
    for case, values in zip(cases, rows):
        if len(values) != 3 or values[0] != case["name"]:
            raise ValueError("порядок результатов не совпадает с набором тестов")
        for value in values[1:]:
            number = float(value)
            if not math.isfinite(number) or (number < 0 and number != -1):
                raise ValueError(f"некорректное время для {case['name']}")
        if (float(values[1]) == -1) != (float(values[2]) == -1):
            raise ValueError("при отказе оба времени должны быть -1")
        row.extend(values[1:])
    results = Path(results)
    results.parent.mkdir(parents=True, exist_ok=True)
    header = result_header(cases)
    descriptor = os.open(results, os.O_RDWR | os.O_CREAT, 0o666)
    with os.fdopen(descriptor, "r+", encoding="utf-8", newline="") as target:
        fcntl.flock(target, fcntl.LOCK_EX)
        exists, legacy, history, contents = _read_open_results(target, header)
        buffer = io.StringIO(newline="")
        writer = csv.writer(buffer, lineterminator="\n")
        if not exists or legacy:
            writer.writerow(header)
        if legacy:
            for previous in history:
                writer.writerow([previous[0], *("unknown" for _key in SYSTEM_COLUMNS),
                                 *previous[1:]])
        writer.writerow(row)
        if legacy:
            backup = results.with_name(results.name + ".before-system.bak")
            version = 1
            while True:
                try:
                    with backup.open("x", encoding="utf-8", newline="") as saved:
                        saved.write(contents)
                        saved.flush()
                        os.fsync(saved.fileno())
                    break
                except FileExistsError:
                    backup = results.with_name(results.name + f".before-system.{version}.bak")
                    version += 1
            target.seek(0)
        else:
            target.seek(0, os.SEEK_END)
        target.write(buffer.getvalue())
        if legacy:
            target.truncate()
        target.flush()
        os.fsync(target.fileno())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    system = modes.add_parser("system")
    system.add_argument("--output", type=Path, required=True)
    run = modes.add_parser("run")
    run.add_argument("--log", type=Path, required=True)
    run.add_argument("--expect", choices=sorted(EXPECTATIONS), required=True)
    for argument in ("n", "m", "s"):
        run.add_argument(f"--{argument}", type=int, required=True)
    run.add_argument("command", nargs=argparse.REMAINDER)
    for mode in ("check", "append"):
        subparser = modes.add_parser(mode)
        subparser.add_argument("--results", type=Path, required=True)
        subparser.add_argument("--manifest", type=Path, required=True)
        if mode == "append":
            subparser.add_argument("--measurements", type=Path, required=True)
            subparser.add_argument("--timestamp", required=True)
            subparser.add_argument("--system-info", type=Path)
    args = parser.parse_args()
    try:
        if args.mode == "system":
            info = collect_system_info()
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(info, ensure_ascii=False, indent=2) + "\n",
                                   encoding="utf-8")
            print(f"Система: {info['system']}; архитектура: {info['architecture']}; "
                  f"CPU: {info['cpu']}; компьютер: {info['hostname']}")
        elif args.mode == "run":
            command = args.command[1:] if args.command[:1] == ["--"] else args.command
            if not command:
                raise ValueError("не задана команда запуска")
            t1, t2, status, passed = execute_case(command, args.log,
                                                 args.expect, args.n, args.m, args.s)
            print(f"{t1}\t{t2}\t{status}\t{int(passed)}")
        else:
            cases = read_manifest(args.manifest)
            if args.mode == "check":
                check_results(args.results, cases)
            else:
                info = (json.loads(args.system_info.read_text(encoding="utf-8"))
                        if args.system_info is not None else None)
                append_results(args.results, cases, args.measurements, args.timestamp, info)
    except (OSError, ValueError, csv.Error) as error:
        parser.exit(2, f"Ошибка: {error}\n")


if __name__ == "__main__":
    main()
