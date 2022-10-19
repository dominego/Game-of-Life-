import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import math

# setting up the values for the grid
ON = 255
OFF = 0


class Grid:
    def __init__(
        self,
        size=100,
        pGame=1,
        denGeneration=0.5,
        density=0,
        grid=np.array([]),
        vals=[OFF, ON],
        upBord=[],
        downBord=[],
        leftBord=[],
        rightBord=[],
        frameNum=500,
    ):
        if size>8:
            self.size = size
        else:
            self.size=8
        self.pGame = pGame
        self.denGeneration = denGeneration
        self.density = density
        self.grid = np.zeros((size, size))
        self.fig, self.ax = plt.subplots()
        self.img = self.ax.imshow(self.grid, interpolation="nearest")
        self.vals = vals
        self.upBord = upBord
        self.downBord = downBord
        self.leftBord = leftBord
        self.rightBord = rightBord
        self.frameNum = frameNum

    def pbc(self):
        for i in range(self.size):
            self.upBord.append(self.grid[i][self.size - 1])
            self.downBord.append(self.grid[i][0])
            self.leftBord.append(self.grid[self.size - 1][i])
            self.rightBord.append(self.grid[0][i])

    def offbc(self):
        for i in range(self.size):
            self.upBord.append(self.vals[0])
            self.downBord.append(self.vals[0])
            self.leftBord.append(self.vals[0])
            self.rightBord.append(self.vals[0])

    def onbc(self):
        for i in range(self.size):
            self.upBord.append(self.vals[1])
            self.downBord.append(self.vals[1])
            self.leftBord.append(self.vals[1])
            self.rightBord.append(self.vals[1])

    # Returns a grid of NxN random values with defined density "den"
    def randomGrid(self):
        self.grid = np.random.choice(
            self.vals,
            self.size * self.size,
            p=[1 - self.denGeneration, self.denGeneration],
        ).reshape(self.size, self.size)

    def gridWithBc(self, BoundaryConditions):
        M = self.size + 2
        gridTot = np.random.choice(self.vals, M * M, p=[1, 0]).reshape(M, M)
        self.randomGrid()
        grid = self.grid
        gridTot[1 : 1 + self.size, 1 : 1 + self.size] = grid
        BoundaryConditions()
        gridTot[0, 1 : self.size + 1] = self.upBord
        gridTot[M - 1, 1 : self.size + 1] = self.downBord
        gridTot[1 : self.size + 1, 0] = self.leftBord
        gridTot[1 : self.size + 1, M - 1] = self.rightBord
        if BoundaryConditions() == self.onbc():
            gridTot[0, 0] = ON
            gridTot[M - 1, M - 1] = ON
            gridTot[0, M - 1] = ON
            gridTot[M - 1, 0] = ON
        self.grid = gridTot

    # Returns the same grid with a random generated MxM subgrid with top-left corner in (i,j)
    # if the subgrid should be put ouside the grid, the initial grid is not changed
    def addSubGrid(self, i, j, M):
        if (i+M) < self.size and (j+M) < self.size and i>0 and j>0:
            subcell = np.random.choice(
                self.vals, M * M, p=[self.denGeneration, 1 - self.denGeneration]
            ).reshape(M, M)
            self.grid[i : i + M, j : j + M] = subcell

    # Add a glider in position (i,j), if i and j are given such that the glider should be put outside of the grid
    # it is not put
    def addGlider(self, i, j):
        if (i+3) < self.size and (j+3) < self.size and i>0 and j>0:
            glider = np.array([[0, 0, 255], [255, 0, 255], [0, 255, 255]])
            self.grid[i : i + 3, j : j + 3] = glider

    # Add a Gosper glider gun in postion (i,j), if i and j are given such that the gosper glider gun should
    # be put outside of the grid it is not put
    def addGosperGliderGun(self, i, j):
        if (i+11) < self.size and (j+38) < self.size and i>0 and j>0:
            gun = np.zeros(11 * 38).reshape(11, 38)

            gun[5][1] = gun[5][2] = 255
            gun[6][1] = gun[6][2] = 255

            gun[3][13] = gun[3][14] = 255
            gun[4][12] = gun[4][16] = 255
            gun[5][11] = gun[5][17] = 255
            gun[6][11] = gun[6][15] = gun[6][17] = gun[6][18] = 255
            gun[7][11] = gun[7][17] = 255
            gun[8][12] = gun[8][16] = 255
            gun[9][13] = gun[9][14] = 255

            gun[1][25] = 255
            gun[2][23] = gun[2][25] = 255
            gun[3][21] = gun[3][22] = 255
            gun[4][21] = gun[4][22] = 255
            gun[5][21] = gun[5][22] = 255
            gun[6][23] = gun[6][25] = 255
            gun[7][25] = 255

            gun[3][35] = gun[3][36] = 255
            gun[4][35] = gun[4][36] = 255

            self.grid[i : i + 11, j : j + 38] = gun

    # Returns the density of the grid
    def density(self):
        alive = 0
        death = 0
        for i in range(self.size):
            for j in range(self.size):
                if self.grid[i, j] == ON:
                    alive += 1
                else:
                    death += 1
            density = alive / (alive + death)
        self.density = density

    # Divide the grid in rectangular grid lx x ly and returns an array with the density of each regions
    # it check if lx and ly divide the grid exactly, if not it returns an empty array
    def densityRegions(self, lx, ly):
        densityRegions = []
        if self.size % lx == 0 and self.size % ly == 0 and lx>0 and ly>0:
            for n in range(int(self.size / lx)):
                for m in range(int(self.size / ly)):
                    alive = 0
                    death = 0
                    for i in range(0 + n * lx, lx + n * lx):
                        for j in range(0 + m * ly, ly + m * ly):
                            if self.grid[i, j] == ON:
                                alive += 1
                            else:
                                death += 1
                    densityRegions.append(alive / (alive + death))

        return densityRegions

    # Return the density of a regions lx x ly with top-left corner in i,j
    # if the region selected goes outside the grid it returns -1
    def densityZone(self, i, j, lx, ly):
        alive = 0
        death = 0
        if (i + lx) < self.size and (j + ly) < self.size and lx>0 and ly>0:
            for i in range(i, i + lx):
                for j in range(j, j + ly):
                    if self.grid[i, j] == ON:
                        alive += 1
                    else:
                        death += 1
            densityZone = alive / (alive + death)

        else:
            densityZone = -1

        return densityZone

    # Returns three densities: density1,density0,densityUncoupled
    ## density 1 is the number of couples of near alive cells along x and y divided by the number of cells
    ## density 0 is the number of couples of near death cells along x and y divided by the number of cells
    ## densityUncopled is the number of couples of near cells in different states along x and y divided by the number of cells
    ### Note: this method dose not look all the possibile couples along x or y because it shifts by two cells at the time.
    ###       it dose not look diagonaly couple. For N big enough this approximation can be done for statistical reason
    def densityCouple(self):
        densityZonesx = []
        densityZonesy = []
        for n in range(int(self.size / 2)):
            for m in range(int(self.size / 1)):
                alive = 0
                death = 0
                for i in range(0 + n * 2, 2 + n * 2):
                    for j in range(0 + m * 1, 1 + m * 1):
                        if self.grid[i, j] == ON:
                            alive += 1
                        else:
                            death += 1
                densityZonesx.append(alive / (alive + death))
        couples1x = filter(lambda couple: couple == 1, densityZonesx)
        couples0x = filter(lambda couple: couple == 0, densityZonesx)
        uncouplesx = filter(lambda couple: couple == 0.5, densityZonesx)

        for n in range(int(self.size / 1)):
            for m in range(int(self.size / 2)):
                alive = 0
                death = 0
                for i in range(0 + n * 1, 1 + n * 1):
                    for j in range(0 + m * 2, 2 + m * 2):
                        if self.grid[i, j] == ON:
                            alive += 1
                        else:
                            death += 1
                densityZonesy.append(alive / (alive + death))
        couples1y = filter(lambda couple: couple == 1, densityZonesy)
        couples0y = filter(lambda couple: couple == 0, densityZonesy)
        uncouplesy = filter(lambda couple: couple == 0.5, densityZonesy)

        density1 = ((len(list(couples1x)) + len(list(couples1y)))) / (
            self.size * self.size
        )
        density0 = (len(list(couples0x)) + len(list(couples0y))) / (
            self.size * self.size
        )
        densityUncoupled = (len(list(uncouplesx)) + len(list(uncouplesy))) / (
            self.size * self.size
        )
        return density1, density0, densityUncoupled

    # Update the state of the grid taking care of the periodic boundary conditions
    def updateBC(self, frameNum, inf=0, sup=0, balance=0):
        # copy grid since the update happen simultaneously
        # for calculation and go line by line
        newGrid = self.grid.copy()
        for i in range(1, self.size + 1):
            for j in range(1, self.size + 1):

                # compute 8-neighbor sum
                # using toroidal boundary conditions - x and y wrap around
                # so that the simulaton takes place on a toroidal surface.
                # %N is done because when i,j is on the edge of the grid, some of their neighbours are in the boundary
                # region.
                total = int(
                    (
                        self.grid[i % self.size, (j - 1) % self.size]
                        + self.grid[i % self.size, (j + 1) % self.size]
                        + self.grid[(i - 1) % self.size, j % self.size]
                        + self.grid[(i + 1) % self.size, j % self.size]
                        + self.grid[(i - 1) % self.size, (j - 1) % self.size]
                        + self.grid[(i - 1) % self.size, (j + 1) % self.size]
                        + self.grid[(i + 1) % self.size, (j - 1) % self.size]
                        + self.grid[(i + 1) % self.size, (j + 1) % self.size]
                    )
                    / 255
                )

                # apply Conway's rules with probability pGame and an alternative rule with probability 1-pGame
                if random.uniform(0, 1) > float(1 - self.pGame):
                    if self.grid[i, j] == ON:
                        if (total < 2) or (total > 3):

                            newGrid[i, j] = OFF
                    else:
                        if total == 3:
                            newGrid[i, j] = ON
                else:
                    self.altRules(newGrid, i, j, total, inf, sup, balance)

        # update data
        self.grid[:] = newGrid[:]
        self.img = self.ax.imshow(self.grid, interpolation="nearest")

    # Update the state of the grid without PBC
    def update(self, frameNum, inf=0, sup=0, balance=0):
        # copy grid since the update happen simultaneously
        # for calculation and we go line by line
        newGrid = self.grid.copy()
        for i in range(1, self.size - 1):
            for j in range(1, self.size - 1):

                total = int(
                    (
                        self.grid[i, (j - 1)]
                        + self.grid[i, (j + 1)]
                        + self.grid[(i - 1), j]
                        + self.grid[(i + 1), j]
                        + self.grid[(i - 1), (j - 1)]
                        + self.grid[(i - 1), (j + 1)]
                        + self.grid[(i + 1), (j - 1)]
                        + self.grid[(i + 1), (j + 1)]
                    )
                    / 255
                )

                # apply Conway's rules apply Conway's rules with probability pGame and an alternative rule with probability 1-pGame
                if random.uniform(0, 1) > float(1 - self.pGame):
                    if self.grid[i, j] == ON:
                        if (total < 2) or (total > 3):

                            newGrid[i, j] = OFF
                    else:
                        if total == 3:
                            newGrid[i, j] = ON
                else:
                    self.altRules(newGrid, i, j, total, inf, sup, balance)
        # update data
        self.grid[:] = newGrid[:]
        self.img = self.ax.imshow(self.grid, interpolation="nearest")

    #The altRule update the grid in (i,j) using a Conway's-like rule but the theshold value are imposed with 
    # inf, sup and balance
    def altRules(self, newGrid, i, j, total, inf=0, sup=0, balance=0):
        if inf == sup == balance == 0:
            newGrid[i, j] = random.choice(self.vals)
        else:
            if self.grid[i, j] == ON:
                if (total < inf) or (total > sup):

                    newGrid[i, j] = OFF
            else:
                if total == balance:
                    newGrid[i, j] = ON
        return newGrid

    #This perform an animation of "frameNum" instants and initialize the grid with randomGrid() so the initial density 
    # will be almost equal to denGeneration. The value inf,sup and baalnce can be passed to define an alternative
    # rule of the game
    def animation(self, inf=0, sup=0, balance=0):
        self.randomGrid()
        ani = animation.FuncAnimation(
            self.fig,
            self.update,
            frames=self.frameNum,
            fargs=(inf, sup, balance),
            interval=0.0001,
            save_count=1,
            repeat=False,
        )
        plt.show()
    #This perform an animation of "frameNum" instants taking care of the boundary conditions. The grid is initialised 
    # using gridWithBc(BoundaryConditions), thus it will have the boundary conditions selected from onbc(),offbc() and
    # pbc() while the rest of the grid is generated with randomGrid()
    def animationBC(self, BoundaryConditions, inf=0, sup=0, balance=0):
        self.gridWithBc(BoundaryConditions)
        ani = animation.FuncAnimation(
            self.fig,
            self.updateBC,
            frames=self.frameNum,
            fargs=(inf, sup, balance),
            interval=0.0001,
            save_count=1,
            repeat=False,
        )
        plt.show()

