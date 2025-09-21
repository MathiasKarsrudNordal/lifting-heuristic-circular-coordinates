### **Lifting Cocycles: From Heuristic to Theory**

--- 
**Repository Structure**

```bash
├─ figs       # Generated figures
├─ notebooks  # Example of Algorithm 1 in practice (Example 4.3) and code to reproduce Figure 3
├─ scripts    # Scripts to reproduce results and plots
└─ src        # Reusable functions and core logic
```

**Get Started** 

Clone the repository and activate the Julia environment.

```Julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Scripts in `scripts/` are provided to regenerate the figures from the paper. Figures are written to the `figs/` directory by default.
