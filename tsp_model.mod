
# TSP Model using Miller-Tucker-Zemlin (MTZ) Constraints

set CITIES;  # Set of cities
param DISTANCE{CITIES, CITIES} >= 0;  # Distance matrix

# Decision variables
var x{CITIES, CITIES} binary;  # 1 if the tour goes from city i to j, 0 otherwise
var u{CITIES} >= 0;  # Subtour elimination variables

# Objective: Minimize total travel distance
minimize Total_Distance:
    sum {i in CITIES, j in CITIES: i != j} DISTANCE[i, j] * x[i, j];

# Each city must be visited exactly once
subject to VisitOnce{i in CITIES}:
    sum {j in CITIES: i != j} x[i, j] = 1;

# Each city must be departed from exactly once
subject to DepartOnce{j in CITIES}:
    sum {i in CITIES: i != j} x[i, j] = 1;

# Subtour elimination constraints
subject to SubtourElimination{i in CITIES, j in CITIES: i != j}:
    u[i] - u[j] + card(CITIES) * x[i, j] <= card(CITIES) - 1;

# Prevent self-loops
subject to NoSelfLoop{i in CITIES}:
    x[i, i] = 0;
