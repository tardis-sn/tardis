(function(){
  'use strict';
  /**
   * Tag input directive for use in list fields
   * based on http://codepen.io/webmatze/pen/isuHh
   */
  angular
      .module('tardisApp')
      .directive('file', file);

  function file(){
    var directive = {
      scope: {
        file: '='
      },
      link: link
    }
    return directive; 

    function link(scope, elem, attrs){
      elem.bind('change', function(event){
        var files = event.target.files;
        var file = files[0];
        scope.file = file ? file.name : undefined;
        scope.$apply();
      });
    }
  }

})();