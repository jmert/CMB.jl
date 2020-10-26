export PolarizationConventions, StokesCrossFields

module PolarizationConventions
    export IAUConv, HealpixConv

    """
        @enum Convention IAUConv HealpixConv

    An enumeration to specify the two types of polarization conventions used to describe
    Stokes Q/U coordinate systems.
    """
    @enum Convention IAUConv HealpixConv
end

module StokesCrossFields
    using BitFlags
    export TT, QT, UT, TQ, QQ, UQ, TU, QU, UU, NO_FIELD, TPol, Pol

    """
        @bitflag Field TT QT UT TQ QQ UQ TU QU UU NO_FIELD=0

    A bitfield for identifying combinations of Stokes fields. There are 9 subblocks, named
    as the Cartesian product of elements T, Q, and U:

        TT  TQ  TU
        QT  QQ  QU
        UT  UQ  UU
    """
    @bitflag Field TT QT UT TQ QQ UQ TU QU UU NO_FIELD=0

    @inline function count_pos(b, i::Int)
        isset = b & (one(b) << (i - 1)) != zero(b)
        psum = count_ones(b & ((one(b) << i) - 1))
        return ifelse(isset, psum, zero(psum))
    end
    field_offsets(f::Field) = ntuple(i -> count_pos(Unsigned(f), i), 9)
    field_count(f::Field) = count_ones(Unsigned(f))

    function minrow(fields::Field)
        F = Unsigned(fields)
        f = (F | (F >> 0x3) | (F >> 0x6)) & 0x07
        return trailing_zeros(f) + 0x1
    end
    function maxrow(fields::Field)
        F = Unsigned(fields)
        f = (F | (F >> 0x3) | (F >> 0x6)) & 0x07
        return 8sizeof(F) - leading_zeros(f)
    end

    function mincol(fields::Field)
        F = Unsigned(fields)
        f = (F | (F >> 0x1) | (F >> 0x2)) & 0b001001001 #= 0x49 =#
        f = (f | (f >> 0x2) | (f >> 0x4)) & 0b111 #= 0x07 =#
        return trailing_zeros(f) + 0x1
    end
    function maxcol(fields::Field)
        F = Unsigned(fields)
        f = (F | (F >> 0x1) | (F >> 0x2)) & 0b001001001 #= 0x49 =#
        f = (f | (f >> 0x2) | (f >> 0x4)) & 0b111 #= 0x07 =#
        return 8sizeof(F) - leading_zeros(f)
    end

    """
        const TPol = QT | UT | TQ | TU

    An alias for the temperature-cross-polarization Stokes field combinations.
    """
    const TPol = QT | UT | TQ | TU

    """
        const Pol  = QQ | UQ | QU | UU

    An alias for the polarization-only sub-blocks Stokes field combinations.
    """
    const Pol  = QQ | UQ | QU | UU
end

