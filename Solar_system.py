from skyfield.api import load, Topos
from datetime import datetime, timedelta
# from datetime import timedelta
import numpy as np
from skyfield.almanac import find_discrete, risings_and_settings
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Coordinates import *

def convert_to_hms(value):
    """Converts a decimal hour or degree value to hours, minutes, seconds."""
    hours = int(value)
    minutes = int((value - hours) * 60)
    seconds = (value - hours - minutes/60) * 3600
    return hours, abs(minutes), abs(seconds)

def get_celestial_positions(time):
    # Load ephemeris data
    eph = load('de421.bsp')
    earth, sun, moon = eph['earth'], eph['sun'], eph['moon']

    # Calculate the position of the Sun relative to Earth at the given time
    astrometric_sun = earth.at(time).observe(sun)
    apparent_sun = astrometric_sun.apparent()
    ra_sun, dec_sun, _ = apparent_sun.radec()

    # Calculate the position of the Moon relative to Earth at the given time
    astrometric_moon = earth.at(time).observe(moon)
    apparent_moon = astrometric_moon.apparent()
    ra_moon, dec_moon, _ = apparent_moon.radec()

    return {
        "sun": (ra_sun.hours- 1/60, dec_sun.degrees +4/60),
        "moon": (ra_moon.hours -8/60, dec_moon.degrees)
    }

def conv_eq_to_cart(ra, dec, distance):
    x = np.cos(np.deg2rad(ra))*np.cos(np.deg2rad(dec))*distance
    y = np.sin(np.deg2rad(ra))*np.cos(np.deg2rad(dec))*distance
    z = np.sin(np.deg2rad(dec))*distance
    return [x, y, z]

def find_direction_cosines_plane(points):

    p1, p2, p3 = points
    # Convert points to numpy arrays for vector operations
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    # Compute two vectors in the plane
    v1 = p2 - p1
    v2 = p3 - p1

    # Compute the normal vector using the cross product
    normal_vector = np.cross(v1, v2)

    # Normalize the normal vector to get the direction cosines
    magnitude = np.linalg.norm(normal_vector)
    direction_cosines = normal_vector / magnitude

    return direction_cosines

def main():

    gl = [0, 22, 120]
    gb = [0, 0, 0]
    galactic_ra, galactic_dec = conv_gal_to_eq(gl, gb)
    print(galactic_ra, galactic_dec)
    sv = conv_eq_to_cart(galactic_ra, galactic_dec, 1)
    x, y, z = sv
    print(list(zip(x, y, z)))

    dc = find_direction_cosines_plane(list(zip(x, y, z)))
    a, b, c = dc
    d =0
    print(dc)


    # simulation starts from current time to one full orbit
    start = np.datetime64(datetime.now())#+ timedelta( days = -7, hours =+16))
    # new_time = start_time 
    
    # Convert numpy.datetime64 back to datetime for compatibility with Skyfield
    start_time = start.astype('M8[ms]').astype(datetime)

    # Calculate Sun's and Moon's positions at the start time
    ts = load.timescale()
    print(start, start_time, ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second))

    positions = get_celestial_positions(ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second))
    
    sun_ra, sun_dec = positions["sun"]
    moon_ra, moon_dec = positions["moon"]

    solar_sv = conv_eq_to_cart(sun_ra*15, sun_dec, 1)
    lunar_sv = conv_eq_to_cart(moon_ra*15, moon_dec, 1)

# Convert RA and Dec to hours, minutes, seconds
    sun_ra_hms = convert_to_hms(sun_ra)
    sun_dec_dms = convert_to_hms(sun_dec)
    
    moon_ra_hms = convert_to_hms(moon_ra)
    moon_dec_dms = convert_to_hms(moon_dec)

    # Print in hours, minutes, seconds format

    print(f"\n Sun RA, Dec: {sun_ra*15}° {sun_dec}°")
    print(fr"\n Moon Ra Dec: {moon_ra*15}° {moon_dec}°")

    print(f"Sun Position - RA: {sun_ra_hms[0]}h {sun_ra_hms[1]}m {sun_ra_hms[2]:.2f}s, Dec: {sun_dec_dms[0]}° {sun_dec_dms[1]}' {sun_dec_dms[2]:.2f}\"")
    print(f"Moon Position - RA: {moon_ra_hms[0]}h {moon_ra_hms[1]}m {moon_ra_hms[2]:.2f}s, Dec: {moon_dec_dms[0]}° {moon_dec_dms[1]}' {moon_dec_dms[2]:.2f}\"")

    print(f"\n Sun Direction Cosines: {solar_sv}")
    print(f"\n Moon Direction Cosines: {lunar_sv}")
    
    # Create 3D scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    distance = 100
    # Plot the vector (arrow) using quiver
    # Scatter plot
    earth = ax.scatter(0, 0, 0, color='blue', s = 0.6*distance, label = 'Earth')
    sun = ax.quiver(solar_sv[0]*distance, solar_sv[1]*distance, solar_sv[2]*distance, -solar_sv[0], -solar_sv[1], -solar_sv[2], length=0.3*distance, color='orange', arrow_length_ratio=0.5, label= 'Sun light')
    moon = ax.quiver(lunar_sv[0]*distance, lunar_sv[1]*distance, lunar_sv[2]*distance, -lunar_sv[0], -lunar_sv[1], -lunar_sv[2], length=0.3*distance, color='grey', arrow_length_ratio=0.5, label = 'Moon light')
    
    print("Earth--------",type(earth))
    print("Sun--------",type(sun))
    # Set labels for the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # ax.set_zscale("log")

    xlim = ylim = zlim = (-distance, distance)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)

    num_points = 1000

    # Calculate and plot intercepts on the XY plane (z=0)
    if c != 0:
        X = np.linspace(xlim[0], xlim[1], num_points)
        Y = (d - a*X +  distance*c) / b
        mask = (Y <=  distance) & (Y >= - distance)
        X = X[mask] 
        Y = Y[mask]     
        Z = np.ones_like(X)*- distance
        ax.plot(X, Y, Z, color='purple' )

        X = np.linspace(xlim[0], xlim[1], num_points)
        Y = (d - a*X -  distance*c) / b
        mask = (Y <=  distance) & (Y >= - distance)
        X = X[mask] 
        Y = Y[mask]     
        Z = np.ones_like(X)* distance
        ax.plot(X, Y, Z, color='purple' )
        
    if b != 0:
        Z = np.linspace(zlim[0], zlim[1], num_points)
        X = (d - c*Z  -  distance*b) / a
        mask = (X <=  distance) & (X >= - distance)
        Z = Z[mask]
        X = X[mask]
        Y = np.ones_like(X)* distance
        ax.plot(X, Y, Z, color='purple' )

        Z = np.linspace(zlim[0], zlim[1], num_points)
        X = (d - c*Z  +  distance*b) / a
        mask = (X <=  distance) & (X >= - distance)
        Z = Z[mask]
        X = X[mask]
        Y = np.ones_like(X)*- distance
        ax.plot(X, Y, Z, color='purple' )

    if a != 0:
        Y =  np.linspace(ylim[0], ylim[1], num_points)
        Z = (d - b*Y +  distance*a) / c
        mask = (Z <=  distance) & (Z >= - distance)

        Z = Z[mask] 
        Y = Y[mask]
        X = np.ones_like(Y)*- distance
        ax.plot(X, Y, Z, color='purple', label = 'Galactic plane' )

        Y =  np.linspace(ylim[0], ylim[1], num_points)
        Z = (d - b*Y -  distance*a) / c
        mask = (Z >= - distance) & (Z<= distance)
        Z = Z[mask] 
        Y = Y[mask]     
        X = np.ones_like(Y)* distance
        ax.plot(X, Y, Z, color='purple' )

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Show the plot

    plt.legend()
    plt.show()

    
    return

if __name__ == '__main__':
    main()
