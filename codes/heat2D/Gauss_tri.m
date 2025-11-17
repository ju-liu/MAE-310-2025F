function [qp, qw] = Gauss_tri(num_pts)
%GAUSS_TRIANGLE  Gauss quadrature on the reference triangle.
%
%   [qp, qw] = gauss_triangle(num_pts)
%
%   Input:
%       num_pts : number of quadrature points
%                 3, 4, 6, 13, 19, 37
%
%   Output:
%       qp : [num_pts x 3]  (r,s,t),
%            r + s + t = 1
%       qw : [num_pts x 1] weights

qp = zeros(num_pts, 3);
qw = zeros(num_pts, 1);

switch num_pts
    case 3
        qp(1,:) = [0.5, 0.5, 0.0];
        qp(2,:) = [0.5, 0.0, 0.5];
        qp(3,:) = [0.0, 0.5, 0.5];
        qw(:)   = 1.0 / 3.0;

    case 4
        qp(1,:) = [1/3, 1/3, 1/3];
        qw(1)   = -0.5625;

        w = 0.520833333333333;
        qp(2,:) = [0.6, 0.2, 0.2];
        qp(3,:) = [0.2, 0.6, 0.2];
        qp(4,:) = [0.2, 0.2, 0.6];
        qw(2:4) = w;

    case 6
        a = 0.816847572980459;
        b = 0.091576213509771;
        w = 0.109951743655322;

        qp(1,:) = [a, b, b];
        qp(2,:) = [b, a, b];
        qp(3,:) = [b, b, a];
        qw(1:3) = w;

        a = 0.108103018168070;
        b = 0.445948490915965;
        w = 0.223381589678011;

        qp(4,:) = [a, b, b];
        qp(5,:) = [b, a, b];
        qp(6,:) = [b, b, a];
        qw(4:6) = w;

    case 13
        % 1-point (center)
        qp(1,:) = [1/3, 1/3, 1/3];
        qw(1)   = -0.149570044467670;

        % 3-point symmetric set
        a = 0.479308067841923;
        b = 0.260345966079038;
        w = 0.175615257433204;
        qp(2,:) = [a, b, b];
        qp(3,:) = [b, a, b];
        qp(4,:) = [b, b, a];
        qw(2:4) = w;

        % 3-point symmetric set
        a = 0.869739794195568;
        b = 0.065130102902216;
        w = 0.053347235608839;
        qp(5,:) = [a, b, b];
        qp(6,:) = [b, a, b];
        qp(7,:) = [b, b, a];
        qw(5:7) = w;

        % 6-point general set
        a = 0.638444188569809;
        b = 0.312865496004875;
        c = 1.0 - a - b;
        w = 0.077113760890257;
        qp(8,:)  = [a, b, c];
        qp(9,:)  = [a, c, b];
        qp(10,:) = [b, a, c];
        qp(11,:) = [c, a, b];
        qp(12,:) = [b, c, a];
        qp(13,:) = [c, b, a];
        qw(8:13) = w;

    case 19
        qp(1,:) = [1/3, 1/3, 1/3];
        qw(1)   = 0.037861091200315;

        % 3-point symmetric set
        a = 0.797426985353087;
        b = 0.101286507323456;
        w = 0.037620425413183;
        qp(2,:) = [a, b, b];
        qp(3,:) = [b, a, b];
        qp(4,:) = [b, b, a];
        qw(2:4) = w;

        % 3-point symmetric set
        a = 0.059715871789770;
        b = 0.470142064105115;
        w = 0.078357352244117;
        qp(5,:) = [a, b, b];
        qp(6,:) = [b, a, b];
        qp(7,:) = [b, b, a];
        qw(5:7) = w;

        % 3-point symmetric set
        a = 0.535795346449899;
        b = 0.232102326775050;
        w = 0.116271479656966;
        qp(8,:)  = [a, b, b];
        qp(9,:)  = [b, a, b];
        qp(10,:) = [b, b, a];
        qw(8:10) = w;

        % 3-point symmetric set
        a = 0.941038278231121;
        b = 0.029480860884440;
        w = 0.013444267375166;
        qp(11,:) = [a, b, b];
        qp(12,:) = [b, a, b];
        qp(13,:) = [b, b, a];
        qw(11:13) = w;

        % 6-point general set
        a = 0.738416812340510;
        b = 0.232102326775050;
        c = 0.029480860884440;
        w = 0.037509722455232;
        qp(14,:) = [a, b, c];
        qp(15,:) = [a, c, b];
        qp(16,:) = [b, a, c];
        qp(17,:) = [c, a, b];
        qp(18,:) = [b, c, a];
        qp(19,:) = [c, b, a];
        qw(14:19) = w;

    case 37
        idx = 1;
        qp(idx,:) = [1/3, 1/3, 1/3];
        qw(idx)   = 0.051739766065744133555179145422;

        % 3-point symmetric set
        a = 0.950275662924105565450352089520;
        b = 0.024862168537947217274823955239;
        w = 0.008007799555564801597804123460;
        for k = 1:3
            idx = idx + 1;
            if k == 1, qp(idx,:) = [a, b, b];
            elseif k == 2, qp(idx,:) = [b, a, b];
            else,          qp(idx,:) = [b, b, a];
            end
            qw(idx) = w;
        end

        a = 0.171614914923835347556304795551;
        b = 0.414192542538082326221847602214;
        w = 0.046868898981821644823226732071;
        for k = 1:3
            idx = idx + 1;
            if k == 1, qp(idx,:) = [a, b, b];
            elseif k == 2, qp(idx,:) = [b, a, b];
            else,          qp(idx,:) = [b, b, a];
            end
            qw(idx) = w;
        end

        a = 0.539412243677190440263092985511;
        b = 0.230293878161404779868453507244;
        w = 0.046590940183976487960361770070;
        for k = 1:3
            idx = idx + 1;
            if k == 1, qp(idx,:) = [a, b, b];
            elseif k == 2, qp(idx,:) = [b, a, b];
            else,          qp(idx,:) = [b, b, a];
            end
            qw(idx) = w;
        end

        a = 0.772160036676532561750285570113;
        b = 0.113919981661733719124857214943;
        w = 0.031016943313796381407646220131;
        for k = 1:3
            idx = idx + 1;
            if k == 1, qp(idx,:) = [a, b, b];
            elseif k == 2, qp(idx,:) = [b, a, b];
            else,          qp(idx,:) = [b, b, a];
            end
            qw(idx) = w;
        end

        a = 0.009085399949835353883572964740;
        b = 0.495457300025082323058213517632;
        w = 0.010791612736631273623178240136;
        for k = 1:3
            idx = idx + 1;
            if k == 1, qp(idx,:) = [a, b, b];
            elseif k == 2, qp(idx,:) = [b, a, b];
            else,          qp(idx,:) = [b, b, a];
            end
            qw(idx) = w;
        end

        a = 0.062277290305886993497083640527;
        b = 0.468861354847056503251458179727;
        w = 0.032195534242431618819414482205;
        for k = 1:3
            idx = idx + 1;
            if k == 1, qp(idx,:) = [a, b, b];
            elseif k == 2, qp(idx,:) = [b, a, b];
            else,          qp(idx,:) = [b, b, a];
            end
            qw(idx) = w;
        end

        a = 0.851306504174348550389457672223;
        b = 0.022076289653624405142446876931;
        c = 0.126617206172027096933163647918;
        w = 0.015445834210701583817692900053;
        permsABC = [a b c;
                    a c b;
                    b a c;
                    b c a;
                    c a b;
                    c b a];
        for k = 1:6
            idx       = idx + 1;
            qp(idx,:) = permsABC(k,:);
            qw(idx)   = w;
        end

        a = 0.018620522802520968955913511549;
        b = 0.689441970728591295496647976487;
        c = 0.291937506468887771754472382212;
        w = 0.017822989923178661888748319485;
        permsABC = [a b c;
                    a c b;
                    b a c;
                    b c a;
                    c a b;
                    c b a];
        for k = 1:6
            idx       = idx + 1;
            qp(idx,:) = permsABC(k,:);
            qw(idx)   = w;
        end

        a = 0.267625659273967961282458816185;
        b = 0.635867859433872768286976979827;
        c = 0.096506481292159228736516560903;
        w = 0.037038683681384627918546472190;
        permsABC = [a b c;
                    a c b;
                    b a c;
                    b c a;
                    c a b;
                    c b a];
        for k = 1:6
            idx       = idx + 1;
            qp(idx,:) = permsABC(k,:);
            qw(idx)   = w;
        end

    otherwise
        error('gauss_triangle: num_pts = %d not implemented.', num_pts);
end

qw = 0.5 * qw;
end