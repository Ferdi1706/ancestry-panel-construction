This code provides an optimization algorithm to select optimal markers for determination of ancestry. We use the delta values of the markers as the main quality metric of the markers and aim to find a set that distinguishes all continents evenly good.
Therefore we frame the problem as a combinatorial optimization problem and use the intlinprog() function, a mixed-integer linear programming (MILP) solver.
At the end, we obtain a marker set where the markers are selected in such a way that the sum of the delta values for each continent combination is maximized evenly, thereby ensuring that the selection provides a good separation of the continents.
