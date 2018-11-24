# from csunpos import sunpos
# import pyximport
#
# pyximport.install()

from sunpos.csunpos import sunpos

def main(args):
    az, zen, ra, dec, h = sunpos(args.t, args.lat, args.lon, args.elev, args.temp, args.p, args.dt, args.rad)
    if args.csv:
        #machine readable
        print('{t}, {dt}, {lat}, {lon}, {elev}, {temp}, {p}, {az}, {zen}, {ra}, {dec}, {h}'.format(t=args.t, dt=args.dt, lat=args.lat, lon=args.lon, elev=args.elev,temp=args.temp, p=args.p,az=az, zen=zen, ra=ra, dec=dec, h=h))
    else:
        dr='deg'
        if args.rad:
            dr='rad'
        print("Computing sun position at T = {t} + {dt} s".format(t=args.t, dt=args.dt))
        print("Lat, Lon, Elev = {lat} deg, {lon} deg, {elev} m".format(lat=args.lat, lon=args.lon, elev=args.elev))
        print("T, P = {temp} C, {press} mbar".format(temp=args.temp, press=args.p))
        print("Results:")
        print("Azimuth, zenith = {az} {dr}, {zen} {dr}".format(az=az,zen=zen,dr=dr))
        print("RA, dec, H = {ra} {dr}, {dec} {dr}, {h} {dr}".format(ra=ra, dec=dec, h=h, dr=dr))

if __name__ == '__main__':
    from argparse import ArgumentParser
    import datetime, sys
    parser = ArgumentParser(prog='sunposition',description='Compute sun position parameters given the time and location')
    parser.add_argument('--version',action='version',version='%(prog)s 1.0')
    parser.add_argument('--citation',dest='cite',action='store_true',help='Print citation information')
    parser.add_argument('-t,--time',dest='t',type=str,default='now',help='"now" or date and time (UTC) in "YYYY-MM-DD hh:mm:ss.ssssss" format or a (UTC) POSIX timestamp')
    parser.add_argument('-lat,--latitude',dest='lat',type=float,default=51.48,help='latitude, in decimal degrees, positive for north')
    parser.add_argument('-lon,--longitude',dest='lon',type=float,default=0.0,help='longitude, in decimal degrees, positive for east')
    parser.add_argument('-e,--elevation',dest='elev',type=float,default=0,help='elevation, in meters')
    parser.add_argument('-T,--temperature',dest='temp',type=float,default=14.6,help='temperature, in degrees celcius')
    parser.add_argument('-p,--pressure',dest='p',type=float,default=1013.0,help='atmospheric pressure, in millibar')
    parser.add_argument('-dt',type=float,default=0.0,help='difference between earth\'s rotation time (TT) and universal time (UT1)')
    parser.add_argument('-r,--radians',dest='rad',action='store_true',help='Output in radians instead of degrees')
    parser.add_argument('--csv',dest='csv',action='store_true',help='Comma separated values (time,dt,lat,lon,elev,temp,pressure,az,zen,RA,dec,H)')
    args = parser.parse_args()
    if args.cite:
        print("Implementation: Samuel Bear Powell, 2016")
        print("Algorithm:")
        print("Ibrahim Reda, Afshin Andreas, \"Solar position algorithm for solar radiation applications\", Solar Energy, Volume 76, Issue 5, 2004, Pages 577-589, ISSN 0038-092X, doi:10.1016/j.solener.2003.12.003")
        sys.exit(0)
    if args.t == "now":
        args.t = datetime.datetime.utcnow()
    elif ":" in args.t and "-" in args.t:
        try:
            args.t = datetime.datetime.strptime(args.t,'%Y-%m-%d %H:%M:%S.%f') #with microseconds
        except:
            try:
                args.t = datetime.datetime.strptime(args.t,'%Y-%m-%d %H:%M:%S.') #without microseconds
            except:
                args.t = datetime.datetime.strptime(args.t,'%Y-%m-%d %H:%M:%S')
    else:
        args.t = datetime.datetime.utcfromtimestamp(int(args.t))
    main(args)