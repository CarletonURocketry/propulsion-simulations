classdef testing_app < matlab.apps.AppBase
    
    properties (Access = public)
        HRAP_v2022_06_04                matlab.ui.Figure
        AppTitle                        matlab.ui.control.Label
    end


    % Component initialization
    methods (Access = private)
        % Create UIFigure and components
        function createComponents(app)

            % Create HRAP_v2022_06_04 and hide until all components are created
            app.HRAP_v2022_06_04 = uifigure('Visible', 'off');
            app.HRAP_v2022_06_04.Position = [100 100 850 889];
            app.HRAP_v2022_06_04.Name = 'MATLAB App';
            app.HRAP_v2022_06_04.Resize = 'off';

            % Create AppTitle
            app.AppTitle = uilabel(app.HRAP_v2022_06_04);
            app.AppTitle.HorizontalAlignment = 'center';
            app.AppTitle.FontName = 'Times New Roman';
            app.AppTitle.FontSize = 24;
            app.AppTitle.FontWeight = 'bold';
            app.AppTitle.Position = [250 850 353 40];
            app.AppTitle.Text = 'Hybrid Rocket Analysis Program';

            % Show the figure after all components are created
            app.HRAP_v2022_06_04.Visible = 'on';
        end
    end
    
    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = testing_app

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.HRAP_v2022_06_04)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.HRAP_v2022_06_04)
        end
    end

end