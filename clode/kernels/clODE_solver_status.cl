#ifndef CLODE_SOLVER_STATUS_H_
#define CLODE_SOLVER_STATUS_H_

// Shared solver-status codes written by the kernel entrypoints and surfaced by
// the runtime fetch path.
#define SOLVE_STATUS_COMPLETED 0
#define SOLVE_STATUS_MAX_STEPS_REACHED 1
#define SOLVE_STATUS_TERMINAL_EVENT_REACHED 2
#define SOLVE_STATUS_OUTPUT_CAPACITY_REACHED 3
#define SOLVE_STATUS_STEPPER_FAILED -1

static inline int finalizeSolveStatus(
    const int solveStatus,
    const bool reachedSolveEnd,
    const bool exhaustedStepBudget,
    const bool hitTerminalEvent,
    const bool exhaustedOutputCapacity
)
{
    if (solveStatus != SOLVE_STATUS_COMPLETED || reachedSolveEnd)
        return solveStatus;

    if (hitTerminalEvent)
        return SOLVE_STATUS_TERMINAL_EVENT_REACHED;

    if (exhaustedStepBudget)
        return SOLVE_STATUS_MAX_STEPS_REACHED;

    if (exhaustedOutputCapacity)
        return SOLVE_STATUS_OUTPUT_CAPACITY_REACHED;

    return solveStatus;
}

#endif // CLODE_SOLVER_STATUS_H_