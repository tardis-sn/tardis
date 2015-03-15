(function() {
  'use strict';

  angular
      .module('tardisApp')
      .directive('quantityField', quantityField);

     function quantityField (){
      var directive = {
        restrict: 'E',
        templateUrl: '/static/app/input-data/quantity.field.tmpl.html',
        require: 'ngModel',
        replace: true,
        scope: {
          'name': '@',
          'units': '@',
          'help': '@'
        },
        link: link,
        controller: quantityFieldCtrl,
        controllerAs: 'vm',
        bindToController: true

      }
     
      return directive;

      function link (scope, elem, attrs, ngModelCtrl){

        scope.vm.unitArray = unitArray;
        scope.vm.unit = unitArray()[0];


        scope.$watch('vm.unit + vm.value', function() {
          var value = scope.vm.value || 0;
          ngModelCtrl.$setViewValue(value +' '+ scope.vm.unit);
        });
        
        function unitArray() {
          if (scope.vm.units === undefined) {
            return [];
          }
          return scope.vm.units.split(',').filter(function(unit) {
            return unit !== "";
          });
        }

      }
     
    }

    quantityFieldCtrl.$inject = ['$scope'];

    function quantityFieldCtrl($scope){
      var vm = this;
    }


})();