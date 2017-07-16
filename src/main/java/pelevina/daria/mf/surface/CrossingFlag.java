package pelevina.daria.mf.surface;


import java.util.EnumSet;

public enum CrossingFlag {
    NoCrossing, ClockWise, CounterClockWise, SphereCrossing, GoesToInfty, BottomCrossing, WallCrossing;

    public static EnumSet CONTINUE = EnumSet.of(ClockWise, CounterClockWise, SphereCrossing, BottomCrossing);
}
