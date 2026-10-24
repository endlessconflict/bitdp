/* Retired user-space instructions and cycles of a command and its children (Linux perf_event_open).
 * usage: perfcount CMD [ARG...]   prints "instructions=N cycles=M" on stderr after CMD exits;
 * the exit status is CMD's. Counting starts at exec, so the count covers CMD only, from its first
 * instruction to its exit (input reading and output included). */
#define _GNU_SOURCE
#include <linux/perf_event.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/syscall.h>
#include <sys/wait.h>
#include <unistd.h>

static int counter(pid_t pid, uint64_t config) {
    struct perf_event_attr a;
    memset(&a, 0, sizeof a);
    a.size = sizeof a;
    a.type = PERF_TYPE_HARDWARE;
    a.config = config;
    a.disabled = 1;
    a.enable_on_exec = 1;
    a.inherit = 1;
    a.exclude_kernel = 1;
    a.exclude_hv = 1;
    return (int)syscall(SYS_perf_event_open, &a, pid, -1, -1, 0);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: perfcount CMD [ARG...]\n");
        return 2;
    }
    int go[2];
    if (pipe(go) != 0) return 2;
    pid_t pid = fork();
    if (pid == 0) { /* wait until the counters are attached, then exec */
        char c;
        close(go[1]);
        if (read(go[0], &c, 1) != 1) _exit(127);
        execvp(argv[1], argv + 1);
        _exit(127);
    }
    int fi = counter(pid, PERF_COUNT_HW_INSTRUCTIONS), fc = counter(pid, PERF_COUNT_HW_CPU_CYCLES);
    if (fi < 0 || fc < 0) {
        perror("perf_event_open");
        return 2;
    }
    close(go[0]);
    if (write(go[1], "x", 1) != 1) return 2;
    close(go[1]);
    int st;
    waitpid(pid, &st, 0);
    uint64_t ni = 0, nc = 0;
    if (read(fi, &ni, sizeof ni) != sizeof ni || read(fc, &nc, sizeof nc) != sizeof nc) return 2;
    fprintf(stderr, "instructions=%llu cycles=%llu\n", (unsigned long long)ni, (unsigned long long)nc);
    return WIFEXITED(st) ? WEXITSTATUS(st) : 2;
}
