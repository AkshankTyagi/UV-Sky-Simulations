(xlim[0], xlim[1], num_points)
        Y = (d - a*X +  distance*c) / b
        mask = (Y <=  distance) & (Y >= - distance)
        X = X[mask] 
        Y = Y[mask]     
        Z = np.ones_like(X)*- distance
        ax.plot(X, Y, Z, color='purple' )